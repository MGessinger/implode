#include "padic_ode_solver.h"

int main ()
{
	int return_value = EXIT_SUCCESS;
	slong p, prec, degree, order;

	padic_ctx_t ctx;
	padic_ode_t ODE;
	padic_t num;
	padic_poly_t result;
	flint_rand_t state;

	/* Initialization */
	flint_randinit(state);

	for (slong iter = 0; iter < 100; iter++)
	{
		p = n_randprime(state, 8, 1);
		prec = 2 + n_randint(state, 62);
		degree = 2 + n_randint(state, 8);
		order = n_randint(state, degree);
		if (order <= 0)
			order = 2;

		padic_ctx_init(ctx, &p, 0, prec, PADIC_SERIES);
		padic_init2(num, prec);
		padic_poly_init2(result, 32, prec);

		padic_ode_init_blank(ODE, degree, order, prec);

		/* Setup */
		for (slong i = 0; i <= order(ODE); i++)
		{
			for (slong j = i+1; j <= degree(ODE); j++)
				padic_randtest(padic_ode_coeff(ODE, i, j), state, ctx);
		}
		padic_one(padic_ode_coeff(ODE, order(ODE), order(ODE)));
		padic_set_si(num, order(ODE) - 1, ctx);

		padic_poly_one(result);
		_padic_ode_solve_frobenius(result, ODE, num, 32, ctx);
		padic_poly_shift_left(result, result, order(ODE) - 1, ctx);

		int solved = padic_ode_solves(ODE, result, 30, prec, ctx);
		if (!solved)
			return_value = EXIT_FAILURE;

		padic_clear(num);
		padic_ode_clear(ODE);
		padic_poly_clear(result);
		padic_ctx_clear(ctx);
	}
	flint_randclear(state);
	flint_cleanup();
	return return_value;
}
