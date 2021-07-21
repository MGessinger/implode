#ifndef SOLVER_FUCHS_H_
#define SOLVER_FUCHS_H_

#include "padic_ode.h"

slong padic_ode_valuation (padic_ode_t ODE);
slong padic_ode_solve_fuchs (padic_poly_t res, padic_ode_t ODE, slong num_of_coeffs, padic_ctx_t ctx, slong prec);

#endif /* SOLVER_FUCHS_H_ */
