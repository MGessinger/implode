#include "padic_ode_solver.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong p, prec, n, degree, order, s_rho;

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
		n = 2 + n_randint(state, 62);
		prec = 30 + n_randint(state, 128);
		degree = 3 + n_randint(state, 7);
		order = 2 + n_randint(state, degree-2);

		padic_ctx_init(ctx, &p, 0, prec, PADIC_SERIES);
		padic_init2(rho, prec);

		s_rho = order - 2;
		padic_set_si(rho, s_rho, ctx);

		padic_ode_init_blank(ODE, degree, order, prec);
		padic_ode_solution_init(sol, rho, 2, 0, ctx);

		/* Setup */
		for (slong i = 0; i <= order(ODE); i++)
			for (slong j = i+1; j <= degree(ODE); j++)
				padic_randtest(padic_ode_coeff(ODE, i, j), state, ctx);
		padic_one(padic_ode_coeff(ODE, order, order));
		padic_one(padic_ode_coeff(ODE, order-1, order-1));

		padic_ode_solve_frobenius(sol, ODE, n, ctx);
		padic_poly_shift_left(sol->gens, sol->gens, s_rho, ctx);
		int solved = padic_ode_solves(ODE, sol->gens + 0, n, prec/2, ctx);

		padic_ode_clear(ODE);
		padic_ode_solution_clear(sol);
		padic_clear(rho);
		padic_ctx_clear(ctx);

		if (!solved)
		{
			return_value = EXIT_FAILURE;
			break;
		}
	}
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
