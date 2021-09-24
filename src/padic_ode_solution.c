#include "padic_ode.h"

void padic_ode_solution_init (padic_ode_solution_t sol, padic_t rho, slong mul, slong alpha, padic_ctx_t ctx)
{
	sol->mul = mul;
	sol->M = mul + alpha;

	padic_init2(sol->rho, padic_get_prec(rho));
	padic_set(sol->rho, rho, ctx);

	sol->gens = flint_malloc(sol->M * sizeof(padic_poly_struct));
	for (slong i = 0; i < sol->M; i++)
		padic_poly_init2(sol->gens + i, 16, padic_get_prec(rho));
}

void padic_ode_solution_clear (padic_ode_solution_t sol)
{
	padic_clear(sol->rho);
	for (slong i = 0; i < sol->M; i++)
		padic_poly_clear(sol->gens + i);
	flint_free(sol->gens);
}

void padic_ode_solution_dump (padic_ode_solution_t sol, padic_ctx_t ctx)
{
	flint_printf("Solution adjoint to the exponent "); padic_print(sol->rho, ctx); flint_printf(" of multiplicity %w.\n", sol->mul);
	for (slong i = 0; i < sol->M; i++)
	{
		flint_printf("log(x)^%w *\t", sol->M - 1 - i);
		padic_poly_print_pretty(sol->gens + i, "x", ctx);
		flint_printf("\n\n");
	}
	flint_printf("\n");
}

void _padic_ode_solution_update (padic_ode_solution_t sol, padic_poly_t f, padic_ctx_t ctx)
{
	padic_struct *F;
	padic_t temp1, temp2;

	F = flint_malloc(sol->M * sizeof(padic_struct));
	if (F == NULL)
		return;

	padic_init(temp1);
	padic_init(temp2);

	padic_one(temp2);

	for (slong k = 0; k < sol->M; k++)
	{
		padic_init2(F + k, padic_poly_prec(f));
		padic_poly_evaluate_padic(F + k, f, sol->rho, ctx);
		padic_poly_derivative(f, f, ctx);

		padic_mul(F + k, F + k, temp2, ctx);
		padic_set_si(temp1, sol->M - 1 - k, ctx);
		padic_mul(temp2, temp2, temp1, ctx);
		padic_set_si(temp1, k + 1, ctx);
		padic_div(temp2, temp2, temp1, ctx);
	}

	for (slong n = sol->M - 1; n >= 0; n--)
	{
		padic_poly_scalar_mul_padic(sol->gens + n, sol->gens + n, F + 0, ctx);
		padic_set_si(temp2, n, ctx);
		for (slong k = 1; k <= n; k++)
		{
			padic_poly_scalar_mul_padic(f, sol->gens + (n - k), F + k, ctx);
			padic_poly_add(sol->gens + n, sol->gens + n, f, ctx);

			padic_set_si(temp1, n - k, ctx);
			padic_mul(F + k, F + k, temp1, ctx);
			padic_div(F + k, F + k, temp2, ctx);
		}
		padic_clear(F + n);
	}

	padic_poly_zero(f);

	padic_clear(temp1);
	padic_clear(temp2);
	flint_free(F);
}

void _padic_ode_solution_extend (padic_ode_solution_t sol, slong nu, padic_poly_t g_nu, padic_ctx_t ctx)
{
	padic_t temp;
	padic_poly_t der;

	padic_init2(temp, padic_get_prec(sol->rho));
	padic_poly_init2(der, padic_poly_length(g_nu), padic_poly_prec(g_nu));

	padic_poly_set(der, g_nu, ctx);
	for (slong i = 0; i < sol->M; i++)
	{
		padic_poly_evaluate_padic(temp, der, sol->rho, ctx);
		padic_poly_set_coeff_padic(sol->gens + i, nu, temp, ctx);
		padic_poly_derivative(der, der, ctx);
	}
	padic_poly_zero(g_nu);
	padic_poly_clear(der);
	padic_clear(temp);
}

void _padic_ode_solution_normalize (padic_ode_solution_t sol, padic_ctx_t ctx)
{
	padic_t t;
	padic_init2(t, padic_get_prec(sol->rho));

	slong alpha = sol->M - sol->mul;
	padic_poly_get_coeff_padic(t, sol->gens + alpha, 0, ctx);
	padic_inv(t, t, ctx);

	for (slong i = 0; i < sol->M; i++)
		padic_poly_scalar_mul_padic(sol->gens + i, sol->gens + i, t, ctx);

	padic_clear(t);
}
