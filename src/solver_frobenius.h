#ifndef SOLVER_FROBENIUS_H_
#define SOLVER_FROBENIUS_H_

#include "padic_ode.h"

/* =============== Solutions =============== */

typedef struct {
	padic_t rho;		/* root of the indicial equation */
	slong multiplicity;	/* counts the number of roots equal to rho */
	slong alpha;		/* alpha distinguishes roots within a group */
	padic_poly_struct *gens;/* contains the series g_ν and its derivatives */
} padic_ode_solution_struct;

typedef padic_ode_solution_struct padic_ode_solution_t[1];

void	padic_ode_solution_init (padic_ode_solution_t sol, padic_t rho, slong mul, slong alpha, padic_ctx_t ctx);
void	padic_ode_solution_clear (padic_ode_solution_t sol);
void	padic_ode_solution_dump (padic_ode_solution_t sol, padic_ctx_t ctx);

void	padic_ode_solution_evaluate (padic_t res, padic_ode_solution_t sol, padic_t x, slong mu, padic_ctx_t ctx);

void	_padic_ode_solution_extend (padic_ode_solution_t sol, slong nu, padic_poly_t g_nu, padic_ctx_t ctx);
void	_padic_ode_solution_update (padic_ode_solution_t sol, padic_poly_t f, padic_ctx_t ctx);

/* =============== Solvers =============== */

int	indicial_polynomial (padic_poly_t result, padic_ode_t ODE, slong nu, slong shift, padic_ctx_t ctx);
int	indicial_polynomial_evaluate (padic_t result, padic_ode_t ODE, padic_t rho, slong shift, slong nu, padic_ctx_t ctx);

void	_padic_ode_solve_frobenius (padic_poly_t res, padic_ode_t ODE, padic_t rho, slong degree, padic_ctx_t ctx);
void	padic_ode_solve_frobenius (padic_ode_solution_t sol, padic_ode_t ODE, slong degree, padic_ctx_t ctx);

#endif /* SOLVER_FROBENIUS_H_ */
