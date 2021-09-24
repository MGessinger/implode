# IMPLODE - The Implementation of p-adic linear, ordinary differential equations

Welcome to Implode v1.0.
Implode is a library designed to store and solve differential equations using p-adic arithmetic.
To do this, it uses [FLINT](https://flintlib.org)'s `padic_t` data type.
Solutions are computed as (generalized) power series solutions, making use of a recursion relation between the coefficients.

For a similar library using complex arithmetic, see [Cascade](https://github.com/MGessinger/Cascade).

## Installation

This library uses cmake as its build system, so to compile it simply run
```bash
mkdir build && cd build
cmake ..
make -j
```
Then, to install the library on your system, run
```
sudo make install
```

## Examples

The following program computes a power series solution to the differential equation `x^2 y'' + (3x-1)y' + y = 0`.
This equation is solved by a power series with coefficients `a_n = n!`, which only converges over the p-adic numbers.

```C
#include <implode.h>

int main ()
{
	slong p = 13;
	padic_t rho;
	padic_ctx_t ctx;
	padic_ode_t ode;
	padic_ode_solution_t sol;

	/* Setup */
	padic_init(rho);
	padic_ctx_init(ctx, &p, 0, 20, PADIC_TERSE);
	padic_ode_solution_init(sol, rho, 1, 0, ctx);

	padic_ode_init_blank(ode, 2, 2, 20);
	padic_set_si(padic_ode_coeff(ode, 2, 2), 1, ctx);
	padic_set_si(padic_ode_coeff(ode, 1, 1), 3, ctx);
	padic_set_si(padic_ode_coeff(ode, 1, 0), -1, ctx);
	padic_set_si(padic_ode_coeff(ode, 0, 0), 1, ctx);

	/* Solving and Output */
	padic_ode_solve_frobenius(sol, ode, 10, ctx);

	padic_poly_print_pretty(sol->gens, "x", ctx);
	flint_printf("\n");

	/* Cleanup */
	padic_ode_solution_clear(sol);
	padic_ode_clear(ode);
	padic_ctx_clear(ctx);
	padic_clear(rho);
	flint_cleanup();
	return 0;
}
```

To compile this into an executable, run
```bash
gcc test.c -lflint -limplode
```

The output should look something like this:
```bash
(3628800)*x^10+(362880)*x^9+(40320)*x^8+(5040)*x^7+(720)*x^6+(120)*x^5+(24)*x^4+(6)*x^3+(2)*x^2+x+(1)
```

## Memory Management

Because Flint can cache some constants internally, it is recommended to call *flint_cleanup()* at the end of your main program.
This will clear Flint's internal cache, and guarentee a clean output of Memory tools like Valgrind.

## Dependencies

Implode uses [FLINT](http://flintlib.org/) to store p-adic numbers and polynmials.
Therefore FLINT has to be installed in order to compile and use Implode.
Furthermore, any program you write with Implode must be linked against libflint.
