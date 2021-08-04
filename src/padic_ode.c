#include "padic_ode.h"

/* Setup and memory management*/

static inline int max_degree (padic_poly_t *polys, slong order)
{
	slong deg, poly_max_degree = 0;
	while (polys[order] == NULL || padic_poly_is_zero(polys[order]))
		order--;
	if (order <= 0)
		return -1;

	for (slong i = 0; i <= order; i++)
	{
		if (polys[i] == NULL)
			continue;
		deg = padic_poly_degree(polys[i]);
		if (poly_max_degree < deg)
			poly_max_degree = deg;
	}
	return deg;
}

void padic_ode_init_blank (padic_ode_t ODE, slong degree, slong order, slong prec)
{
	/* Prepare the Differential equation for later use */
	ODE->order = order;
	ODE->degree = degree;
	ODE->polys = NULL;
	ODE->alloc = 0;

	if (degree < 0 || order <= 0)
		return;

	ODE->alloc = (order + 1) * (degree + 1);
	ODE->polys = flint_malloc(ODE->alloc * sizeof(padic_struct));

	for (slong i = 0; i < ODE->alloc; i++)
		padic_init2(ODE->polys + i, prec);
}

void padic_ode_init (padic_ode_t ODE, padic_poly_t *polys, slong order, slong prec, padic_ctx_t ctx)
{
	/* Create a differential operator defined by *polys* */
	slong degree = max_degree(polys, order);

	padic_ode_init_blank (ODE, degree, order, prec);
	for (slong i = 0; i <= order(ODE); i++)
	{
		if (polys[i] == NULL)
			continue;
		for (slong j = 0; j < padic_poly_length(polys[i]); j++)
		{
			padic_poly_get_coeff_padic(padic_ode_coeff(ODE, i, j), polys[i], j, ctx);
		}
	}
}

void padic_ode_clear (padic_ode_t ODE)
{
	if (ODE->alloc <= 0)
		return;

	for (slong i = 0; i < ODE->alloc; i++)
		padic_clear(ODE->polys + i);

	flint_free(ODE->polys);
}

void padic_ode_set (padic_ode_t ODE_out, padic_ode_t ODE_in, slong prec, padic_ctx_t ctx)
{
	padic_ode_clear(ODE_out);
	padic_ode_init_blank(ODE_out, order(ODE_in), degree(ODE_in), prec);

	for (slong i = 0; i < ODE_in->alloc; i++)
		padic_set(ODE_out->polys + i, ODE_in->polys + i, ctx);
}

void padic_ode_dump (padic_ode_t ODE, char *file, padic_ctx_t ctx)
{
	/* Dumps the ODE to file. If file is NULL, dump to stdout */
	FILE *out = stdout;
	if (file != NULL)
		out = fopen(file, "w");
	if (out == NULL)
		return;
	flint_fprintf(out, "Order: %w\nDegree: %w\n", order(ODE), degree(ODE));
	for (slong i = 0; i <= order(ODE); i++)
	{
		flint_fprintf(out, "padic_ode_poly(ODE, %w) = ", i);
		for (slong j = 0; j <= degree(ODE); j++)
		{
			padic_fprint(out, padic_ode_coeff(ODE, i, j), ctx);
			flint_fprintf(out, "\t");
		}
		flint_fprintf(out, "\n");
	}
	if (out != stdout)
		fclose(out);
}
