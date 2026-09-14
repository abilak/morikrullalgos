# Mori and Krull algorithms

This is a C++20 implementation of the finite tests described in Sections 4 and 5 of
the accompanying paper. It handles:

- Mori and Krull decisions for additive evaluation monoids in a real quadratic field;
- construction of a finite additive generating set in the positive quadratic cases;
- the cyclic Mori classification from the minimal polynomial;
- the cyclic Krull algorithm, including the generating-index search and a Hilbert-basis
  saturation test.

All root-counting, semigroup, lattice, and quadratic-field calculations use exact
arbitrary-precision integer or rational arithmetic.

The only build dependency is the header-only Boost.Multiprecision library.

## Build and test

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Input

Quadratic generators are pairs `a,b`, meaning `a + b*sqrt(D)`. Rational coordinates
use `p/q`. The implementation requires positive squarefree `D` and generators spanning
the quadratic field.

```sh
./build/morikrull quadratic-mori 5 1/2,1/2
./build/morikrull quadratic-krull 2 0,2 0,3
```

Cyclic input is the minimal polynomial, with coefficients in constant-first order.
The program assumes the polynomial supplied really is minimal over Q.

```sh
./build/morikrull cyclic-mori -1,-1,1
./build/morikrull cyclic-krull -2,0,1
./build/morikrull cyclic-krull -5,2,-2,1
```

For a finite monoid, `alpha^j = [v0, ...]` reports the coordinates of the additive
generator in the basis `1, alpha, ..., alpha^(d-1)`. A failed Krull test also prints a
lattice hole witnessing failure of integer divisor closedness.

## Algorithms

For the quadratic Mori test, the program checks exactly that every input is a weak
Perron algebraic integer and that at least one conjugate is negative. It then enumerates
monomials by total degree until the current additive span is closed under multiplication
by every input generator. This is the terminating search from Theorem 4.2.

For cyclic input, a Sturm sequence counts positive roots exactly. The modulus test uses
the characteristic polynomial of the tensor square of the companion matrix: its roots
are all products of two conjugates, so a positive root above `alpha^2` exists exactly
when some conjugate has modulus above `alpha`. Once finite generation passes, powers
are reduced through the minimal-polynomial recurrence. The first power in the additive
span of its predecessors is the generating index.

The Krull tests identify the quotient group with an integer lattice. The cone is covered
by all full-dimensional simplicial subcones. Lattice points in their half-open fundamental
parallelepipeds give a finite generating set for the cone lattice monoid; deleting
decomposable points leaves its Hilbert basis. Membership of each Hilbert-basis point in
the original affine monoid is decided by bounded exact integer enumeration against a
strictly positive dual-cone functional.
