#include "padic_ode_solver.h"

void indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	padic_poly_zero(result);
	if (nu > degree(ODE))
		return;
	slong prec = padic_poly_prec(result) + order(ODE);

	padic_t temp1;
	padic_poly_t horner;

	padic_init2(temp1, prec);
	padic_poly_init2(horner, 2, prec);

	padic_one(temp1);
	padic_poly_set_coeff_padic(horner, 1, temp1, ctx);

	padic_poly_fit_length(result, order(ODE)+1);

	for (slong lambda = order(ODE); lambda >= 0; lambda--)
	{
		padic_set_si(temp1, shift - lambda, ctx);
		padic_poly_set_coeff_padic(horner, 0, temp1, ctx);
		padic_poly_mul(result, result, horner, ctx);

		if (lambda + nu <= degree(ODE))
		{
			padic_poly_get_coeff_padic(temp1, result, 0, ctx);
			padic_add(temp1, temp1, padic_ode_coeff(ODE, lambda, lambda + nu), ctx);
			padic_poly_set_coeff_padic(result, 0, temp1, ctx);
		}
	}

	padic_clear(temp1);
	padic_poly_clear(horner);
	return;
}

void indicial_polynomial_evaluate (padic_t result, padic_ode_t ODE, slong nu, padic_t rho, slong shift, padic_ctx_t ctx)
{
	if (nu > degree(ODE))
	{
		padic_zero(result);
		return;
	}

	slong prec = padic_get_prec(rho) + order(ODE);

	padic_t temp1;
	padic_init2(temp1, prec);
	padic_set_si(temp1, shift - order(ODE), ctx);
	padic_add(temp1, rho, temp1, ctx);

	padic_t const_one;
	padic_init2(const_one, prec);
	padic_one(const_one);

	padic_zero(result);
	for (slong lambda = order(ODE); lambda >= 0; lambda--)
	{
		padic_mul(result, result, temp1, ctx);
		padic_add(temp1, temp1, const_one, ctx);
		if (lambda + nu <= degree(ODE))
			padic_add(result, result, padic_ode_coeff(ODE, lambda, lambda+nu), ctx);
	}
	padic_clear(temp1);
	padic_clear(const_one);
	return;
}

void _padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong sol_degree, padic_ctx_t ctx)
{
	slong prec = padic_get_prec(rho) + sol_degree;
	padic_t g_new, temp, g_i;
	padic_init2(g_new, prec);
	padic_init2(temp, prec);
	padic_init2(g_i, prec);

	padic_poly_fit_length(res, sol_degree+1);
	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		padic_zero(g_new);

		slong i = clamp(degree(ODE), 1, nu);
		indicial_polynomial_evaluate(temp, ODE, i, rho, nu-i, ctx);
		do
		{
			padic_poly_get_coeff_padic(g_i, res, nu-i, ctx);
			padic_mul(temp, temp, g_i, ctx);
			padic_sub(g_new, g_new, temp, ctx);

			i--;
			indicial_polynomial_evaluate(temp, ODE, i, rho, nu-i, ctx);
		} while (i > 0);
		padic_div(g_new, g_new, temp, ctx);
		padic_poly_set_coeff_padic(res, nu, g_new, ctx);
	}
	padic_clear(g_new);
	padic_clear(temp);
	padic_clear(g_i);
}

void padic_ode_solve_frobenius (padic_ode_solution_t sol, padic_ode_t ODE, slong sol_degree, padic_ctx_t ctx)
{
	if (sol->multiplicity == 1 && sol->alpha == 0)
	{
		_padic_ode_solve_frobenius(sol->gens + 0, ODE, sol->rho, sol_degree, ctx);
		return;
	}

	padic_poly_t indicial;
	padic_poly_t g_new;
	padic_poly_struct *g_rho;

	g_rho = flint_malloc( (degree(ODE)+1) * sizeof(padic_poly_struct) );
	if (g_rho == NULL)
		return;

	padic_poly_init(indicial);
	padic_poly_init(g_new);
	for (slong i = 0; i <= degree(ODE); i++)
		padic_poly_init(g_rho + i);

	padic_poly_one(g_new);
	padic_poly_set(g_rho + 0, g_new, ctx);
	_padic_ode_solution_extend(sol, 0, g_new, ctx);

	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		/* Compute the new coefficient (as a function of rho) */
		slong i = clamp(degree(ODE), 1, nu+1);
		indicial_polynomial(indicial, ODE, i, nu-i, ctx);
		do
		{
			padic_poly_mul(indicial, indicial, g_rho + i, ctx);
			padic_poly_sub(g_new, g_new, indicial, ctx);

			i--;
			indicial_polynomial(indicial, ODE, i, nu-i, ctx);
		} while (i > 0);

		/* Multiply all relevant g_nu(rho) by f(rho + nu) */
		for (slong i = 1; i <= degree(ODE); i++)
		{
			padic_poly_mul(g_rho + i, g_rho + (i - 1), indicial, ctx);
		}
		padic_poly_set(g_rho + 0, g_new, ctx);

		/* Update the G^(i) */
		_padic_ode_solution_update(sol, indicial, ctx);

		/* Extend the g^(k)_nu */
		_padic_ode_solution_extend(sol, nu, g_new, ctx);
	}

	padic_poly_clear(g_new);
	padic_poly_clear(indicial);
	for (slong i = 0; i <= degree(ODE); i++)
		padic_poly_clear(g_rho + i);
	flint_free(g_rho);
}
