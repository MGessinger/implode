#ifndef PADIC_ODE_H_
#define PADIC_ODE_H_

#include <flint/padic.h>
#include <flint/padic_poly.h>

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
