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
	padic_set_si(padic_ode_coeff(ODE, 1, 2), 7, ctx);
	padic_set_si(padic_ode_coeff(ODE, 2, 2), 1, ctx);

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
	padic_ode_clear(ODE);
	padic_ctx_clear(ctx);
	flint_cleanup();
	return return_value;
}
