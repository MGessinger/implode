#include "padic_ode.h"
#include "solver_frobenius.h"

int main () {
	slong p = 7;
	padic_ctx_t ctx; padic_ctx_init(ctx, &p, 0, 16, PADIC_SERIES);

	/* Initialize the ODE */
	padic_ode_t ODE;
	padic_ode_init_blank(ODE, 2, 2, 20);

	padic_set_si(padic_ode_coeff(ODE, 0, 0), 1, ctx);
	padic_set_si(padic_ode_coeff(ODE, 1, 1), 2, ctx);
	padic_set_si(padic_ode_coeff(ODE, 2, 2), 1, ctx);

	/* Test the computation of indicial equation */
	padic_poly_t res; padic_poly_init(res);
	indicial_polynomial(res, ODE, 0, 0, ctx);
	int return_value = 0;
	if (padic_poly_degree(res) != 2)
		return_value = 1;

	/* Memory Cleanup */
	padic_poly_clear(res);
	padic_ode_clear(ODE);
	padic_ctx_clear(ctx);
	flint_cleanup();
	return return_value;
}
