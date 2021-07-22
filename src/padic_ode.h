#ifndef PADIC_DIFFEQ_H_
#define PADIC_DIFFEQ_H_

#include <flint/padic.h>
#include <flint/padic_poly.h>

typedef struct {
	slong order;
	slong degree;
	padic_t *polys;
} padic_ode;

typedef padic_ode* padic_ode_t;

#define degree(ODE) ((ODE)->degree)
#define order(ODE) ((ODE)->order)
#define diff_eq_poly(ODE,i) (((ODE)->polys)[(i)*(degree(ODE)+1)])
#define diff_eq_coeff(ODE,i,j) (&diff_eq_poly(ODE,i)[j])

/* Memory Management */
void padic_ode_init (padic_ode_t ODE);
void padic_ode_clear(padic_ode_t ODE);

/* I/O */
void padic_ode_dump (padic_ode_t ODE, char *file, padic_ctx_t ctx);

static inline slong clamp (slong in, slong min, slong max)
{
	if (in < min)
		return min;
	else if (in > max)
		return max;
	else
		return in;
}

#endif /* PADIC_DIFFEQ_H_ */
