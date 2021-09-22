#include "padic_ode.h"
#include "implode.h"

void assert_equal (const char *errMsg, padic_t exp, padic_t real, padic_ctx_t ctx)
{
	padic_sub(real, exp, real, ctx);
	if (padic_is_zero(real))
		return;

	flint_printf("%s failed in precision %w\n", errMsg, padic_prec(exp));
	flint_printf("Expected: "); padic_print(exp, ctx); flint_printf("\n");
	flint_printf("Error: "); padic_print(real, ctx); flint_printf("\n\n");
	flint_abort();
}

int main () {
	slong p, prec, degree, order, rho;

	padic_ctx_t ctx;
	padic_t num, exp;
	padic_ode_t ODE;
	padic_poly_t poly, indicial;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	for (slong iter = 0; iter < 100; iter++)
	{
		p = n_randprime(state, 8, 1);
		rho = n_randint(state, 5);
		prec = 2 + n_randint(state, 62);
		degree = 2 + n_randint(state, 8);
		order = n_randint(state, degree);
		if (order <= 0)
			order = 2;

		padic_init2(num, prec);
		padic_init2(exp, prec);

		padic_ctx_init(ctx, &p, 0, prec, PADIC_SERIES);
		padic_poly_init2(poly, 16, prec);
		padic_poly_init2(indicial, 16, prec);

		/* Setup */
		padic_ode_init_blank(ODE, degree, order, prec);
		for (slong i = 0; i <= order(ODE); i++)
			for (slong j = i; j <= degree(ODE); j++)
				padic_randtest(padic_ode_coeff(ODE, i, j), state, ctx);

		padic_one(exp);
		padic_poly_set_coeff_padic(poly, rho, exp, ctx);

		padic_ode_apply(poly, ODE, poly, prec, ctx);

		for (slong i = 0; i <= degree(ODE)+1; i++)
		{
			padic_poly_get_coeff_padic(exp, poly, rho+i, ctx);

			/* Directly compute the "indicial coefficient" */
			padic_set_si(num, rho, ctx);
			indicial_polynomial_evaluate(num, ODE, i, num, 0, ctx);
			assert_equal("Direct computation f(rho)", exp, num, ctx);

			/* Compare to indicial polynomial */
			padic_set_si(num, rho, ctx);
			indicial_polynomial(indicial, ODE, i, 0, ctx);
			padic_poly_evaluate_padic(num, indicial, num, ctx);
			assert_equal("Polynomial computation f(r)|r=rho", exp, num, ctx);

			/* Compute f(0 + rho) instead of f(rho + 0) */
			padic_zero(num);
			indicial_polynomial_evaluate(num, ODE, i, num, rho, ctx);
			assert_equal("Direct computation f(0+rho)", exp, num, ctx);

			/* Again, compare to indicial polynomial */
			padic_zero(num);
			indicial_polynomial(indicial, ODE, i, rho, ctx);
			padic_poly_evaluate_padic(num, indicial, num, ctx);
			assert_equal("Polynomial computation f(r+rho)|r=0", exp, num, ctx);
		}

		/* Memory Cleanup */
		padic_poly_clear(poly);
		padic_poly_clear(indicial);
		padic_clear(num);
		padic_clear(exp);
		padic_ctx_clear(ctx);
		padic_ode_clear(ODE);
	}
	flint_randclear(state);
	flint_cleanup();
	return EXIT_SUCCESS;
}
