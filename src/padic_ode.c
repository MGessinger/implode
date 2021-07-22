#include "padic_ode.h"

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
		flint_fprintf(out, "diff_eq_poly(ODE, %w) = ", i);
		for (slong j = 0; j <= degree(ODE); j++)
		{
			padic_fprint(out, diff_eq_coeff(ODE, i, j), ctx);
			flint_fprintf(out, "\t");
		}
		flint_fprintf(out, "\n");
	}
	if (out != stdout)
		fclose(out);
	return;
}
