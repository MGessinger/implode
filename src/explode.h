#ifndef PADIC_ODE_SOLVER_H
#define PADIC_ODE_SOLVER_H

#include "padic_ode.h"

/* ============================== Fuchs Solver ============================== */

slong	padic_ode_solve_fuchs (padic_poly_t res, padic_ode_t ODE, slong num_of_coeffs, padic_ctx_t ctx);

/* ============================ Frobenius Solver ============================ */

void	indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx);
void	indicial_polynomial_evaluate (padic_t result, padic_ode_t ODE, slong nu, padic_t rho, slong shift, padic_ctx_t ctx);

void	_padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong degree, padic_ctx_t ctx);
void	padic_ode_solve_frobenius (padic_ode_solution_t sol, padic_ode_t ODE, slong degree, padic_ctx_t ctx);

#endif /* PADIC_ODE_SOLVER_H */
