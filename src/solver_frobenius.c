#include "solver_frobenius.h"

void indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx, slong prec)
{
	/* Compute f_ν(ρ-ς) [for definition of f, see Frobenius Paper, equation 3] */
	padic_poly_zero(result);
	padic_poly_fit_length(result, order(ODE)+1);
	if (nu >= degree(ODE))
		return;

	padic_poly_t tmpPoly;
	padic_poly_init(tmpPoly);

	padic_t tmpScalar;
	padic_init2(tmpScalar, prec);

	for (slong i = order(ODE); i >= 0; i--)
	{
		/* Multiply `result` by ((x+shift)-(λ-i)) */
		padic_set_si(tmpScalar, i-shift, ctx);
		padic_poly_scalar_mul_padic(tmpPoly, result, tmpScalar, ctx);
		padic_poly_shift_left(result, result, 1, ctx);
		padic_poly_sub(result, result, tmpPoly, ctx);

		/* Add the coefficient of p_i */
		padic_poly_get_coeff_padic(tmpScalar, result, 0, ctx);
		padic_add(tmpScalar, diff_eq_coeff(ODE, i, i+nu), tmpScalar, ctx);
		padic_poly_set_coeff_padic(result, 0, tmpScalar, ctx);
	}

	padic_clear(tmpScalar);
	padic_poly_clear(tmpPoly);
}

