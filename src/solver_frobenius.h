#ifndef SOLVER_FROBENIUS_H_
#define SOLVER_FROBENIUS_H_

#include "padic_ode.h"

int indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx);
int indicial_polynomial_evaluate (padic_t result, padic_ode_t ODE, padic_t rho, slong shift, slong nu, padic_ctx_t ctx);

void padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong degree, padic_ctx_t ctx);

#endif /* SOLVER_FROBENIUS_H_ */
