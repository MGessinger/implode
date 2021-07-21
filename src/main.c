#include "solver_fuchs.h"
#include "solver_frobenius.h"

int main () {
	slong p = 13;
	slong prec = 64;

	padic_t *poly = malloc(16*sizeof(padic_t));
	if (poly == NULL)
		return -1;
	padic_ctx_t ctx; padic_ctx_init(ctx, &p, 0, prec, PADIC_SERIES);

	/* Initialize the ODE */
	slong coeffs[16] = {0,0,0,0, 1,0,0,0, 0,0,2,0, 0,0,0,1};
	for (int i = 0; i < 16; i++)
	{
		padic_init2(poly[i], prec);
		padic_set_si(poly[i], coeffs[i], ctx);
	}
	padic_ode ODE_struct;
	padic_ode_t ODE = &ODE_struct;
	order(ODE) = 3;
	degree(ODE) = 3;
	ODE_struct.polys = poly;

	padic_ode_dump(ODE, NULL, ctx);

	/* Initialize the power series */
	padic_poly_t res; padic_poly_init2(res, 8, prec);
	padic_t num; padic_init2(num, prec);

	indicial_polynomial(res, ODE, 0, 0, ctx, prec);

	flint_printf("\n\n\nIndicial Equation:\n");
	padic_poly_print(res, ctx);
	
	padic_one(num);
	padic_poly_set_coeff_padic(res, 0, num, ctx);
	padic_poly_set_coeff_padic(res, 1, num, ctx);

	padic_ode_solve_fuchs(res, ODE, 16, ctx, prec);

	flint_printf("\n\n\nPower series solution:\n");
	padic_poly_print(res, ctx);
	flint_printf("\n");

	/* Clear all Storage */
	padic_clear(num);
	padic_poly_clear(res);
	for (int i = 0; i < 16; i++)
		padic_clear(poly[i]);
	free(poly);
	padic_ctx_clear(ctx);
	flint_cleanup();
	return 0;
}
