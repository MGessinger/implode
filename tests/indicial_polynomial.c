#include "padic_ode.h"
#include "solver_frobenius.h"

int main () {
	/* Memory Setup */
	padic_t *poly = malloc(9*sizeof(padic_t));
	if (poly == NULL)
		return -1;

	slong p = 7;
	padic_ctx_t ctx; padic_ctx_init(ctx, &p, 0, 64, PADIC_SERIES);

	/* Initialize the ODE */
	slong coeffs[9] = {1,0,0, 0,2,0, 0,0,1};
	for (int i = 0; i < 9; i++)
	{
		padic_init(poly[i]);
		padic_set_si(poly[i], coeffs[i], ctx);
	}
	padic_ode ODE_struct;
	padic_ode_t ODE = &ODE_struct;
	order(ODE) = 2;
	degree(ODE) = 2;
	ODE_struct.polys = poly;

	/* Test the computation of indicial equation */
	padic_poly_t res; padic_poly_init(res);
	indicial_polynomial(res, ODE, 0, 0, ctx);
	int return_value = 0;
	if (padic_poly_degree(res) != 2)
		return_value = 1;

	/* Memory Cleanup */
	padic_poly_clear(res);
	for (int i = 0; i < 9; i++)
		padic_clear(poly[i]);
	free(poly);
	padic_ctx_clear(ctx);
	flint_cleanup();
	return return_value;
}
