#include "solver_frobenius.h"

int main ()
{
	/* Init */
	padic_ode_solution_t sol;
	padic_ctx_t ctx;
	padic_poly_t f, g;
	padic_poly_struct *H;
	padic_t rho, val;
	slong p = 13;

	flint_rand_t state;
	flint_randinit(state);

	slong n = n_randint(state, 10);
	slong prec = n_randint(state, 64);

	padic_ctx_init(ctx, &p, 0, 16, PADIC_SERIES);
	padic_init2(rho, prec);
	padic_init2(val, prec);
	padic_poly_init2(f, 16, prec);
	padic_poly_init2(g, 16, prec);

	padic_randtest(rho, state, ctx);
	padic_ode_solution_init(sol, rho, n_randint(state, 10), 0, ctx);

	int return_value = EXIT_SUCCESS;

	/* Setup */
	padic_poly_randtest_not_zero(f, state, 20, ctx);

	for (slong i = 0; i < n; i++) {
		padic_poly_init2(H + i, 16, prec);
		padic_poly_randtest_not_zero(g, state, 20, ctx);
		padic_poly_mul(H + i, f, g, ctx);
		_padic_ode_solution_extend(sol, i, g, ctx);
	}

	_padic_ode_solution_update(sol, f, ctx);

	/* Check */
	for (slong j = 0; j < n; j++)
	{
		for (slong i = 0; i < sol->multiplicity; i++)
		{
			padic_poly_evaluate_padic(rho, H + j, sol->rho, ctx);
			padic_poly_get_coeff_padic(val, sol->gens + i, j, ctx);
			padic_sub(rho, rho, val, ctx);
			slong val = padic_get_val(rho);
			if (val != 0 && val < prec - 5) {
				return_value = EXIT_FAILURE;
				break;
			}
			padic_poly_derivative(H + j, H + j, ctx);
		}
		padic_poly_clear(H + j);
	}

	/* Clear */
	flint_randclear(state);
	padic_clear(rho);
	padic_clear(val);
	padic_poly_clear(f);
	padic_poly_clear(g);
	padic_ctx_clear(ctx);
	padic_ode_solution_clear(sol);
	return return_value;
}
