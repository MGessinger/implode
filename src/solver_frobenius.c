#include "solver_frobenius.h"

int indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	padic_poly_zero(result);
	if (nu > degree(ODE))
		return 0;
	slong imax = clamp(degree(ODE) - nu, 0, order(ODE));
	slong prec = padic_poly_prec(result) + imax;

	padic_poly_t tmpPoly;
	padic_poly_init2(tmpPoly, order(ODE)+1, prec);
	padic_poly_fit_length(result, order(ODE)+1);

	padic_t tmpScalar;
	padic_init2(tmpScalar, prec);

	for (; imax >= 0; imax--)
	{
		/* Multiply `result` by ((x+shift)-(λ-i)) */
		padic_set_si(tmpScalar, imax-shift, ctx);
		padic_poly_scalar_mul_padic(tmpPoly, result, tmpScalar, ctx);
		padic_poly_shift_left(result, result, 1, ctx);
		padic_poly_sub(result, result, tmpPoly, ctx);

		/* Add the coefficient of p_i */
		padic_poly_get_coeff_padic(tmpScalar, result, 0, ctx);
		padic_add(tmpScalar, diff_eq_coeff(ODE, imax, imax+nu), tmpScalar, ctx);
		padic_poly_set_coeff_padic(result, 0, tmpScalar, ctx);
	}

	padic_clear(tmpScalar);
	padic_poly_clear(tmpPoly);
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
		padic_add(result, result, diff_eq_coeff(ODE, imax, imax+nu), ctx);
	}
	padic_clear(temp1);
	padic_clear(const_one);
	return 1;
}

void padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong degree, padic_ctx_t ctx)
{
	slong prec = padic_get_prec(rho) + degree;
	padic_poly_fit_length(res, degree);
	padic_t newCoeff, temp, g_nu;
	padic_init2(newCoeff, prec);
	padic_init2(temp, prec);
	padic_init2(g_nu, prec);

	for (slong nu = order(ODE); nu < degree; nu++)
	{
		padic_zero(newCoeff);

		int step = 1;
		int abort = 0;
		while ( !(abort | (step>nu)) )
		{
			abort = indicial_polynomial_evaluate(temp, ODE, rho, nu-step, step, ctx);
			padic_poly_get_coeff_padic(g_nu, res, nu-step, ctx);
			padic_mul(temp, temp, g_nu, ctx);
			padic_sub(newCoeff, newCoeff, temp, ctx);
			step++;
		}
		indicial_polynomial_evaluate(temp, ODE, rho, nu, 0, ctx);
		padic_div(newCoeff, newCoeff, temp, ctx);
		padic_poly_set_coeff_padic(res, nu, newCoeff, ctx);
	}
	padic_clear(newCoeff);
	padic_clear(temp);
	padic_clear(g_nu);
}
