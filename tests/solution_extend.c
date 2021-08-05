#include "solver_frobenius.h"

int main ()
{
	/* Init */
	padic_ode_solution_t sol;
	padic_ctx_t ctx;
	padic_poly_t f, g;
	padic_t rho, val;
	slong p = 13;

	flint_rand_t state;
	flint_randinit(state);

	slong n = n_randint(state, 10);
	int return_value = EXIT_SUCCESS;

	padic_ctx_init(ctx, &p, 0, 16, PADIC_SERIES);
	padic_init(rho);
	padic_init(val);
	padic_poly_init(f);
	padic_poly_init(g);

	/* Test */
	padic_randtest(rho, state, ctx);
	padic_ode_solution_init(sol, rho, n_randint(state, 5), 0, ctx);
	
	for (slong i = 0; i < n; i++)
	{
		padic_poly_randtest_not_zero(f, state, 20, ctx);
		padic_poly_set(g, f, ctx);
		_padic_ode_solution_extend(sol, i, g, ctx);
		for (slong j = 0; j < sol->multiplicity; j++)
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

	flint_randclear(state);
	flint_cleanup();
	return return_value;
}