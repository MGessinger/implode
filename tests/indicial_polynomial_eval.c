#include "padic_ode.h"
#include "solver_frobenius.h"

int main () {
	int return_value = EXIT_SUCCESS;
	slong p, shift, nu, degree, order;

	padic_ctx_t ctx;
	padic_ode_t ODE;
	padic_t num, exp;
	padic_poly_t poly;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	p = 7;
	order = n_randint(state, 10);
	degree = order + n_randint(state, 5);
	shift = n_randint(state, 10);
	nu = n_randint(state, degree);

	padic_ctx_init(ctx, &p, 0, PADIC_DEFAULT_PREC, PADIC_SERIES);
	padic_init(num);
	padic_init(exp);
	padic_poly_init(poly);

	padic_ode_init_blank(ODE, degree, order, PADIC_DEFAULT_PREC);

	/* Setup */
	for (slong i = 0; i <= order(ODE); i++)
		for (slong j = i; j <= degree(ODE); j++)
			padic_randtest(padic_ode_coeff(ODE, i, j), state, ctx);

	padic_randtest(num, state, ctx);
	indicial_polynomial(poly, ODE, nu, shift, ctx);
	padic_poly_evaluate_padic(exp, poly, num, ctx);

	indicial_polynomial_evaluate(num, ODE, num, shift, nu, ctx);

	flint_printf("num: "); padic_print(num, ctx); flint_printf("\n");
	flint_printf("exp: "); padic_print(exp, ctx); flint_printf("\n");

	if (!padic_equal(num, exp))
		return_value = EXIT_FAILURE;

	/* Memory Cleanup */
	padic_clear(num);
	padic_clear(exp);
	padic_ode_clear(ODE);
	padic_ctx_clear(ctx);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
