#include "padic_ode.h"
#include "solver_frobenius.h"

int main () {
	int return_value = EXIT_SUCCESS;
	slong p, degree, order, rho;

	padic_ctx_t ctx;
	padic_t num, exp;
	padic_ode_t ODE;
	padic_poly_t poly, indicial;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	p = n_randprime(state, 8, 1);
	degree = 1 + n_randint(state, 9);
	order = 2 + n_randint(state, degree-2);
	rho = n_randint(state, 5);

	padic_init(num);
	padic_init(exp);

	padic_ctx_init(ctx, &p, 0, PADIC_DEFAULT_PREC, PADIC_SERIES);
	padic_poly_init(poly);
	padic_poly_init(indicial);

	/* Setup */
	padic_ode_init_blank(ODE, degree, order, PADIC_DEFAULT_PREC);
	for (slong i = 0; i <= order(ODE); i++)
		for (slong j = i; j <= degree(ODE); j++)
			padic_randtest(padic_ode_coeff(ODE, i, j), state, ctx);

	padic_one(exp);
	padic_poly_set_coeff_padic(poly, rho, exp, ctx);

	padic_ode_apply(poly, ODE, poly, PADIC_DEFAULT_PREC * 3 / 2, ctx);

	for (slong i = 0; i <= degree(ODE)+1; i++)
	{
		padic_poly_get_coeff_padic(exp, poly, rho+i, ctx);

		/* Directly compute the "indicial coefficient" */
		padic_set_si(num, rho, ctx);
		indicial_polynomial_evaluate(num, ODE, i, num, 0, ctx);
		if (!padic_equal(num, exp))
		{
			return_value = EXIT_FAILURE;
			break;
		}

		/* Compare to indicial polynomial */
		padic_set_si(num, rho, ctx);
		indicial_polynomial(indicial, ODE, i, 0, ctx);
		padic_poly_evaluate_padic(num, indicial, num, ctx);
		if (!padic_equal(num, exp))
		{
			return_value = EXIT_FAILURE;
			break;
		}

		/* Compute f(0 + rho) instead of f(rho + 0) */
		padic_zero(num);
		indicial_polynomial_evaluate(num, ODE, i, num, rho, ctx);
		if (!padic_equal(num, exp))
		{
			return_value = EXIT_FAILURE;
			break;
		}

		/* Again, compare to indicial polynomial */
		padic_zero(num);
		indicial_polynomial(indicial, ODE, i, rho, ctx);
		padic_poly_evaluate_padic(num, indicial, num, ctx);
		if (!padic_equal(num, exp))
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}

	/* Memory Cleanup */
	padic_poly_clear(poly);
	padic_poly_clear(indicial);
	padic_clear(num);
	padic_clear(exp);
	padic_ctx_clear(ctx);
	padic_ode_clear(ODE);
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
