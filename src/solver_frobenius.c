#include "solver_frobenius.h"

/* =============== Solutions =============== */

void padic_ode_solution_init (padic_ode_solution_t sol, padic_t rho, slong mul, slong alpha, padic_ctx_t ctx)
{
	sol->multiplicity = mul;
	sol->alpha = alpha;

	padic_init2(sol->rho, padic_get_prec(rho));
	padic_set(sol->rho, rho, ctx);

	sol->gens = flint_malloc(mul * sizeof(padic_poly_struct));
	for (slong i = 0; i < mul; i++)
		padic_poly_init2(sol->gens + i, 16, padic_get_prec(rho));
}

void padic_ode_solution_clear (padic_ode_solution_t sol)
{
	padic_clear(sol->rho);
	for (slong i = 0; i < sol->multiplicity; i++)
		padic_poly_clear(sol->gens + i);
	flint_free(sol->gens);
}

void _padic_ode_solution_extend (padic_ode_solution_t sol, slong nu, padic_poly_t g_nu, padic_ctx_t ctx)
{
	padic_t temp;
	padic_init2(temp, padic_get_prec(sol->rho));
	for (slong i = 0; i < sol->multiplicity; i++)
	{
		padic_poly_evaluate_padic(temp, g_nu, sol->rho, ctx);
		padic_poly_set_coeff_padic(sol->gens + i, nu, temp, ctx);
		padic_poly_derivative(g_nu, g_nu, ctx);
	}
	padic_clear(temp);
}

/* =============== Solvers =============== */

int indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	padic_poly_zero(result);
	if (nu > degree(ODE))
		return 0;
	slong imax = clamp(degree(ODE) - nu, 0, order(ODE));
	slong prec = padic_poly_prec(result) + imax;

	padic_poly_fit_length(result, order(ODE)+1);

	padic_t temp1, temp2;
	padic_init2(temp1, prec);
	padic_init2(temp2, prec);

	for (; imax >= 0; imax--)
	{
		for (slong i = padic_poly_length(result); i >= 1; i--)
		{
			padic_set_si(temp1, imax-shift, ctx);
			padic_poly_get_coeff_padic(temp2, result, i, ctx);
			padic_mul(temp2, temp2, temp1, ctx);
			padic_poly_get_coeff_padic(temp1, result, i-1, ctx);
			padic_add(temp2, temp2, temp1, ctx);
			padic_poly_set_coeff_padic(result, i, temp2, ctx);
		}

		/* Add the coefficient of p_i */
		padic_poly_get_coeff_padic(temp2, result, 0, ctx);
		padic_add(temp2, padic_ode_coeff(ODE, imax, imax+nu), temp2, ctx);
		padic_poly_set_coeff_padic(result, 0, temp2, ctx);
	}

	padic_clear(temp1);
	padic_clear(temp2);
	return 1;
}

int indicial_polynomial_evaluate (padic_t result, padic_ode_t ODE, padic_t rho, slong shift, slong nu, padic_ctx_t ctx)
{
	if (nu > degree(ODE))
	{
		padic_zero(result);
		return 0;
	}

	slong imax = clamp(degree(ODE) - nu, 0, order(ODE));
	slong prec = padic_get_prec(rho) + imax;

	padic_t temp1;
	padic_init2(temp1, prec);
	padic_set_si(temp1, imax+shift, ctx);
	padic_sub(temp1, rho, temp1, ctx);

	padic_t const_one;
	padic_init2(const_one, prec);
	padic_one(const_one);

	padic_zero(result);
	for (; imax >= 0; imax--)
	{
		padic_mul(result, result, temp1, ctx);
		padic_add(temp1, temp1, const_one, ctx);
		padic_add(result, result, padic_ode_coeff(ODE, imax, imax+nu), ctx);
	}
	padic_clear(temp1);
	padic_clear(const_one);
	return 1;
}

void _padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong sol_degree, padic_ctx_t ctx)
{
	slong prec = padic_get_prec(rho) + sol_degree;
	padic_t newCoeff, temp, g_nu;
	padic_init2(newCoeff, prec);
	padic_init2(temp, prec);
	padic_init2(g_nu, prec);

	padic_poly_fit_length(res, sol_degree);
	for (slong nu = order(ODE); nu < sol_degree; nu++)
	{
		padic_zero(newCoeff);

		slong shift = 1;
		int abort = 0;
		while ( !(abort | (shift>nu)) )
		{
			abort = indicial_polynomial_evaluate(temp, ODE, rho, nu-shift, shift, ctx);
			padic_poly_get_coeff_padic(g_nu, res, nu-shift, ctx);
			padic_mul(temp, temp, g_nu, ctx);
			padic_sub(newCoeff, newCoeff, temp, ctx);
			shift++;
		}
		indicial_polynomial_evaluate(temp, ODE, rho, nu, 0, ctx);
		padic_div(newCoeff, newCoeff, temp, ctx);
		padic_poly_set_coeff_padic(res, nu, newCoeff, ctx);
	}
	padic_clear(newCoeff);
	padic_clear(temp);
	padic_clear(g_nu);
}

void padic_ode_solve_frobenius (padic_ode_solution_t sol, padic_ode_t ODE, slong sol_degree, padic_ctx_t ctx)
{
	if (sol->multiplicity == 1)
	{
		_padic_ode_solve_frobenius(sol->gens + 0, ODE, sol->rho, sol_degree, ctx);
		return;
	}

	padic_poly_t indicial;
	padic_poly_t g_new;
	padic_poly_struct *g_rho = flint_malloc( (degree(ODE)+1) * sizeof(padic_poly_struct) );
	if (g_rho == NULL)
		return;

	padic_poly_init(indicial);
	padic_poly_init(g_new);
	for (slong i = 0; i <= degree(ODE); i++)
		padic_poly_init(g_rho + i);

	padic_poly_one(g_rho + 0);
	_padic_ode_solution_extend(sol, 0, g_rho + 0, ctx);

	for (slong nu = 1; nu < sol_degree; nu++)
	{
		// pass
	}

	padic_poly_clear(g_new);
	padic_poly_clear(indicial);
	for (slong i = 0; i <= degree(ODE); i++)
		padic_poly_clear(g_rho + i);
	flint_free(g_rho);
}
