#ifndef PADIC_ODE_H_
#define PADIC_ODE_H_

#include <flint/padic.h>
#include <flint/padic_poly.h>

/* ========================= Differential Operators ========================= */

typedef struct {
	slong order;
	slong degree;
	slong alloc;
	padic_struct *polys;
} padic_ode;

typedef padic_ode padic_ode_t[1];

#define degree(ODE) ((ODE)->degree)
#define order(ODE) ((ODE)->order)
#define padic_ode_poly(ODE, i) (((ODE)->polys) + (i)*(degree(ODE)+1))
#define padic_ode_coeff(ODE, i, j) (padic_ode_poly(ODE, i) + (j))

/* Setup and memory management */
void	padic_ode_init_blank (padic_ode_t ODE, slong degree, slong order, slong prec);
void	padic_ode_init (padic_ode_t ODE, padic_poly_t *polys, slong order, slong prec, padic_ctx_t ctx);
void	padic_ode_clear (padic_ode_t ODE);
void	padic_ode_set (padic_ode_t ODE_out, padic_ode_t ODE_in, slong prec, padic_ctx_t ctx);

/* I/O */
void	padic_ode_dump (padic_ode_t ODE, char *file, padic_ctx_t ctx);

/* Transformations */
slong	padic_ode_reduce (padic_ode_t ODE);
slong	padic_ode_valuation (padic_ode_t ODE);

/* Differential Action */
void	padic_ode_apply (padic_poly_t out, padic_ode_t ODE, padic_poly_t in, slong prec, padic_ctx_t ctx);
int	padic_ode_solves (padic_ode_t ODE, padic_poly_t res, slong deg, slong prec, padic_ctx_t ctx);

/* =============================== Solutions ================================ */

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

/* ================================ Inlines ================================= */

static inline slong clamp (slong in, slong min, slong max)
{
	if (in < min)
		return min;
	else if (in > max)
		return max;
	else
		return in;
}

#endif /* PADIC_ODE_H_ */
