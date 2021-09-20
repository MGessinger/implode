#include "padic_ode.h"

int main ()
{
	/* Init */
	int return_value = EXIT_SUCCESS;
	slong p, n, prec;
	padic_ode_solution_t sol;
	padic_ctx_t ctx;
	padic_poly_t f, g;
	padic_t rho, val;

	flint_rand_t state;
	flint_randinit(state);

	for (slong iter = 0; iter < 100; iter++)
	{
		p = n_randprime(state, 8, 1);
		n = n_randint(state, 10);
		prec = 2 + n_randint(state, 62);

		padic_ctx_init(ctx, &p, 0, prec, PADIC_SERIES);
		padic_init2(rho, prec);
		padic_init2(val, prec);
		padic_poly_init2(f, 20, prec);
		padic_poly_init2(g, 20, prec);

		/* Test */
		padic_randtest(rho, state, ctx);
		padic_ode_solution_init(sol, rho, n_randint(state, 5), 0, ctx);

		for (slong i = 0; i < n; i++)
		{
			padic_poly_randtest_not_zero(f, state, 20, ctx);
			padic_poly_set(g, f, ctx);
			_padic_ode_solution_extend(sol, i, g, ctx);
			for (slong j = 0; j < sol->M; j++)
			{
				padic_poly_evaluate_padic(rho, f, sol->rho, ctx);
				padic_poly_get_coeff_padic(val, sol->gens + j, i, ctx);
				if (!padic_equal(rho, val))
					return_value = EXIT_FAILURE;
				padic_poly_derivative(f, f, ctx);
			}
		}

		/* Clear */
		padic_clear(rho);
		padic_clear(val);
		padic_poly_clear(f);
		padic_poly_clear(g);
		padic_ode_solution_clear(sol);
		padic_ctx_clear(ctx);
	}
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
