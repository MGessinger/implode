#include "solver_frobenius.h"

int main ()
{
	padic_ode_solution_t sol;
	padic_ctx_t ctx;
	padic_t rho;
	slong p = 13;

	padic_ctx_init(ctx, &p, 0, 16, PADIC_SERIES);

	padic_init(rho);
	padic_one(rho);

	padic_ode_solution_init(sol, rho, 2, 0, ctx);

	padic_ode_solution_clear(sol);

	padic_ctx_clear(ctx);
	padic_clear(rho);
	return 0;
}
