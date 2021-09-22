#include "padic_ode.h"
#include "implode.h"

void indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	padic_poly_zero(result);
	nu += padic_ode_valuation(ODE);
	if (nu > degree(ODE))
		return;

	slong prec = padic_poly_prec(result) + order(ODE) + 2;

	padic_t temp1;
	padic_poly_t out, horner;

	padic_init2(temp1, prec);
	padic_poly_init2(horner, 2, prec);
	padic_poly_init2(out, order(ODE) + 1, prec);

	padic_one(temp1);
	padic_poly_set_coeff_padic(horner, 1, temp1, ctx);

	slong lambda = clamp(degree(ODE) - nu, 0, order(ODE));
	for (; lambda >= 0; lambda--)
	{
		padic_set_si(temp1, shift - lambda, ctx);
		padic_poly_set_coeff_padic(horner, 0, temp1, ctx);
		padic_poly_mul(out, out, horner, ctx);

		if (lambda + nu < 0)
			continue;

		padic_poly_get_coeff_padic(temp1, out, 0, ctx);
		padic_add(temp1, temp1, padic_ode_coeff(ODE, lambda, lambda + nu), ctx);
		padic_poly_set_coeff_padic(out, 0, temp1, ctx);
	}
	padic_poly_set(result, out, ctx);

	padic_clear(temp1);
	padic_poly_clear(horner);
	padic_poly_clear(out);
}

void indicial_polynomial_evaluate (padic_t result, padic_ode_t ODE, slong nu, padic_t rho, slong shift, padic_ctx_t ctx)
{
	nu += padic_ode_valuation(ODE);
	if (nu > degree(ODE))
	{
		padic_zero(result);
		return;
	}

	slong prec = padic_get_prec(result);

	padic_t temp1, out;
	padic_init2(temp1, prec);
	padic_init2(out, prec);

	slong lambda = clamp(degree(ODE) - nu, 0, order(ODE));
	for (; lambda >= 0; lambda--)
	{
		padic_set_si(temp1, shift - lambda, ctx);
		padic_add(temp1, rho, temp1, ctx);
		padic_mul(out, out, temp1, ctx);

		if (lambda + nu < 0)
			continue;

		padic_add(out, out, padic_ode_coeff(ODE, lambda, lambda + nu), ctx);
	}
	padic_set(result, out, ctx);

	padic_clear(temp1);
	padic_clear(out);
}

void _padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong sol_degree, padic_ctx_t ctx)
{
	slong prec = padic_get_prec(rho) + sol_degree;
	padic_t g_new, indicial, g_i;
	padic_init2(g_new, prec);
	padic_init2(indicial, prec);
	padic_init2(g_i, prec);

	padic_poly_one(res);
	padic_poly_fit_length(res, sol_degree + 1);
	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		padic_zero(g_new);

		slong i = clamp(degree(ODE), 1, nu);
		indicial_polynomial_evaluate(indicial, ODE, i, rho, nu - i, ctx);
		do
		{
			padic_poly_get_coeff_padic(g_i, res, nu - i, ctx);
			padic_mul(indicial, indicial, g_i, ctx);
			padic_sub(g_new, g_new, indicial, ctx);

			i--;
			indicial_polynomial_evaluate(indicial, ODE, i, rho, nu - i, ctx);
		} while (i > 0);
		padic_div(g_new, g_new, indicial, ctx);
		padic_poly_set_coeff_padic(res, nu, g_new, ctx);
	}
	padic_clear(g_new);
	padic_clear(indicial);
	padic_clear(g_i);
}

void padic_ode_solve_frobenius (padic_ode_solution_t sol, padic_ode_t ODE, slong sol_degree, padic_ctx_t ctx)
{
	if (sol->M == 1)
	{
		_padic_ode_solve_frobenius(sol->gens, ODE, sol->rho, sol_degree, ctx);
		return;
	}

	slong prec = padic_prec(sol->rho);

	padic_t temp;
	padic_poly_t indicial;
	padic_poly_t g_new;
	padic_poly_struct *g_rho;

	g_rho = flint_malloc( degree(ODE) * sizeof(padic_poly_struct) );
	if (g_rho == NULL)
		return;

	padic_init2(temp, prec);
	padic_poly_init2(indicial, order(ODE) + 1, prec);
	padic_poly_init2(g_new, order(ODE) + 1, prec);

	for (slong i = 0; i < degree(ODE); i++)
		padic_poly_init2(g_rho + i, order(ODE) + 1, prec);
	padic_poly_one(g_rho);

	for (slong i = 0; i < sol->M; i++)
	{
		padic_poly_zero(sol->gens + i);
		padic_poly_fit_length(sol->gens + i, sol_degree + 1);
	}
	padic_poly_one(sol->gens);

	for (slong nu = 1; nu <= sol_degree; nu++)
	{
		/* Compute the new coefficient (as a function of rho) */
		slong i = clamp(nu, 1, degree(ODE));
		indicial_polynomial(indicial, ODE, i, nu - i, ctx);
		do
		{
			padic_poly_mul(indicial, indicial, g_rho + (i - 1), ctx);
			padic_poly_sub(g_new, g_new, indicial, ctx);

			i--;
			indicial_polynomial(indicial, ODE, i, nu - i, ctx);
		} while (i > 0);

		/* Rescale the indicial polynomial, to keep coefficients small */
		indicial_polynomial_evaluate(temp, ODE, 0, sol->rho, nu, ctx);
		if (!padic_is_zero(temp))
		{
			padic_inv(temp, temp, ctx);
			padic_poly_scalar_mul_padic(indicial, indicial, temp, ctx);
			padic_poly_scalar_mul_padic(g_new, g_new, temp, ctx);
		}

		/* Multiply all relevant g_nu(rho) by f(rho + nu) */
		int all_zero = padic_poly_is_zero(g_new);
		for (i = degree(ODE) - 1; i > 0; i--)
		{
			padic_poly_mul(g_rho + i, g_rho + (i - 1), indicial, ctx);
			all_zero &= padic_poly_is_zero(g_rho + i);
		}
		padic_poly_set(g_rho, g_new, ctx);

		/* Update the G^(i) */
		_padic_ode_solution_update(sol, indicial, ctx);

		/* Extend the g^(k)_nu */
		_padic_ode_solution_extend(sol, nu, g_new, ctx);

		if (all_zero)
			break;
	}

	padic_poly_clear(g_new);
	padic_poly_clear(indicial);
	for (slong i = 0; i < degree(ODE); i++)
		padic_poly_clear(g_rho + i);
	flint_free(g_rho);
	padic_clear(temp);
}
