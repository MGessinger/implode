#include "solver_frobenius.h"

int main ()
{
	slong p = 13;
	slong prec = 128;
	padic_ctx_t ctx;
	padic_ctx_init(ctx, &p, 0, 16, PADIC_SERIES);

	padic_ode ODE_struct;
	padic_ode_t ODE = &ODE_struct;
	order(ODE) = 1;
	degree(ODE) = 1;
	ODE->polys = malloc(4*sizeof(padic_t));
	for (int i = 0; i < 4; i++)
	{
		padic_init2(diff_eq_coeff(ODE, i/2, i%2), prec);
		padic_set_si(diff_eq_coeff(ODE, i/2, i%2), i%2, ctx);
	}

	padic_poly_t res;
	padic_poly_init2(res, 32, prec);
	padic_poly_one(res);

	padic_t rho; padic_init2(rho, prec);
	padic_ode_solve_frobenius(res, ODE, rho, 32, ctx);

	int return_value = 0;
	padic_one(rho);
	padic_t err; padic_init2(err, prec);
	for (int i = 1; i < 32; i++)
	{
		padic_set_si(err, i, ctx);
		padic_div(rho, rho, err, ctx);
		padic_poly_get_coeff_padic(err, res, i, ctx);
		padic_sub(err, rho, err, ctx);

		/* Calculate precision loss:
		 *                 (big)              (small)          (big)	*/
		slong mag_err = (padic_val(err) - padic_val(rho)) - padic_prec(rho);
		return_value |= (mag_err > 1);
	}

	padic_clear(rho);
	padic_clear(err);
	for (int i = 0; i < 4; i++)
		padic_clear(diff_eq_coeff(ODE, i/2,i%2));
	free(ODE->polys);
	padic_ctx_clear(ctx);
	padic_poly_clear(res);
	return return_value;
}
