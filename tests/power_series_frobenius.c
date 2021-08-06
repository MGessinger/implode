#include "solver_frobenius.h"

int main ()
{
	slong p = 13;

	padic_ctx_t ctx;
	padic_t rho;
	padic_ode_t ODE;
	padic_ode_solution_t sol;

	padic_ctx_init(ctx, &p, 0, PADIC_DEFAULT_PREC, PADIC_SERIES);

	padic_init(rho);
	padic_set_si(rho, 2, ctx);
	padic_inv(rho, rho, ctx);
	padic_ode_solution_init(sol, rho, 2, 0, ctx);

	padic_ode_init_blank(ODE, 2, 2, PADIC_DEFAULT_PREC);
	padic_set_si(padic_ode_coeff(ODE, 2, 2), 1, ctx);
	padic_set_si(padic_ode_coeff(ODE, 0, 0), 4, ctx);
	padic_inv(padic_ode_coeff(ODE, 0, 0), padic_ode_coeff(ODE, 0, 0), ctx);

	padic_ode_dump(ODE, NULL, ctx);

	padic_ode_solve_frobenius(sol, ODE, 10, ctx);

	padic_ode_solution_dump(sol, ctx);

	padic_ode_clear(ODE);
	padic_ode_solution_clear(sol);
	padic_clear(rho);
	padic_ctx_clear(ctx);
	return EXIT_SUCCESS;
}
