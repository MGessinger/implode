#include "padic_ode_solver.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong p, prec, degree, order;

	flint_rand_t state;
	padic_ctx_t ctx;
	padic_ode_t ODE;
	padic_t rho;
	padic_ode_solution_t sol;

	/* Initialization */
	flint_randinit(state);

	for (slong iter = 0; iter < 100; iter++)
	{
		p = n_randprime(state, 8, 1);
		prec = 2 + n_randint(state, 30);
		degree = 2 + n_randint(state, 8);
		order = n_randint(state, degree);
		if (order <= 0)
			order = 2;

		padic_ctx_init(ctx, &p, 0, prec, PADIC_SERIES);
		padic_init2(rho, prec);

		padic_ode_init_blank(ODE, degree, order, prec);
		padic_ode_solution_init(sol, rho, 2, 0, ctx);

		/* Setup */
		for (slong i = 0; i <= order(ODE); i++)
			for (slong j = i+1; j <= degree(ODE); j++)
				padic_randtest(padic_ode_coeff(ODE, i, j), state, ctx);

		padic_one(padic_ode_coeff(ODE, order, order));
		padic_set_si(padic_ode_coeff(ODE, order-1, order-1), order-1, ctx);

		padic_ode_solve_frobenius(sol, ODE, 32, ctx);

		int solved = padic_ode_solves(ODE, sol->gens + 0, 32, prec, ctx);
		if (!solved)
			return_value = EXIT_FAILURE;

		padic_ode_clear(ODE);
		padic_ode_solution_clear(sol);
		padic_clear(rho);
		padic_ctx_clear(ctx);
	}
	flint_randclear(state);
	flint_cleanup();
	return EXIT_SUCCESS;
}
