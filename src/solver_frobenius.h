#ifndef SOLVER_FROBENIUS_H_
#define SOLVER_FROBENIUS_H_

#include "padic_ode.h"

void indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx, slong prec);

#endif /* SOLVER_FROBENIUS_H_ */
