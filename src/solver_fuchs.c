#include "solver_fuchs.h"

slong padic_ode_valuation (padic_ode_t ODE)
{
	slong val = 0;
	for (int i = 0; i <= order(ODE); i++)
	{
		slong v = i;
		while (padic_is_zero(diff_eq_coeff(ODE, i, i-v)))
			v--;

		if (v > val)
			val = v;
	}
	return val;
}

slong padic_ode_solve_fuchs (padic_poly_t res, padic_ode_t ODE, slong num_of_coeffs, padic_ctx_t ctx)
{
	/* Iteratively compute the first num_of_coeffs coefficients of the power series solution of the ODE around zero */
	padic_t temp1; padic_init(temp1);
	padic_t temp2; padic_init(temp2);
	padic_t new_coeff; padic_init(new_coeff);

	fmpz_t fac; fmpz_init(fac);
	slong i_min, i_max;
	slong v = padic_ode_valuation(ODE);

	padic_poly_fit_length(res, num_of_coeffs);
	for (slong b_max = v; b_max < num_of_coeffs; b_max++)
	{
		padic_zero(new_coeff);
		slong exp = b_max - v;
		/* Loop through the known coefficients of the power series */
		slong b_min = clamp(exp - degree(ODE), 0, b_max);
		padic_set(temp2, diff_eq_coeff(ODE, 0, exp - b_min), ctx);
		slong b = b_min;
		do {
			padic_poly_get_coeff_padic(temp1, res, b, ctx);
			padic_mul(temp2, temp1, temp2, ctx);
			padic_sub(new_coeff, new_coeff, temp2, ctx);

			padic_zero(temp2);
			fmpz_one(fac);
			/* Loop through the polynomials */
			b++;
			i_min = clamp(b - exp, 0, v);
			i_max = clamp(b - b_min, 0, order(ODE));
			for (slong i = i_min; i <= i_max; i++)
			{
				padic_set_fmpz(temp1, fac, ctx);
				padic_mul(temp1, diff_eq_coeff(ODE, i, i + exp - b), temp1, ctx);
				padic_add(temp2, temp2, temp1, ctx);
				fmpz_mul_si(fac, fac, b-i);
			}
			fmpz_rfac_uiui(fac, b-i_min+1, i_min);
			padic_set_fmpz(temp1, fac, ctx);
			padic_mul(temp2, temp2, temp1, ctx);
		} while (b != b_max);
		padic_div(new_coeff, new_coeff, temp2, ctx);
		padic_poly_set_coeff_padic(res, b_max, new_coeff, ctx);
	}

	padic_clear(new_coeff);
	padic_clear(temp1);
	padic_clear(temp2);
	fmpz_clear(fac);
	return 1;
}
