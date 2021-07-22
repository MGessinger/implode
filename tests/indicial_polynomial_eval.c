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
	slong coeffs[9] = {1,0,0, 0,2,7, 0,0,1};
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

	/* Test the direct evaluation of indicial equation */
	padic_t num; padic_init(num);
	padic_t exp; padic_init(exp);
	padic_set_si(num, 5, ctx);
	padic_set_si(exp, 7*5, ctx);
	indicial_polynomial_evaluate(num, ODE, num, 0, 1, ctx);
	int return_value = 0;
	if (!padic_equal(num, exp))
		return_value = 1;

	/* Memory Cleanup */
	padic_clear(num);
	padic_clear(exp);
	for (int i = 0; i < 9; i++)
		padic_clear(poly[i]);
	free(poly);
	padic_ctx_clear(ctx);
	flint_cleanup();
	return return_value;
}
