# morikrull

Decides Mori and Krull properties for the additive evaluation monoids
`N_0[alpha_1, ..., alpha_n]` of Yau paper: the cyclic case of Section 5 in any degree, and
the quadratic multigenerator case of Section 4.

Needs FLINT, GMP, and Normaliz with its C++ headers.

## Build

```sh
cmake -S . -B build -DNORMALIZ_ROOT=/path/to/normaliz/prefix
cmake --build build
ctest --test-dir build --output-on-failure
```

## Run

```
morikrull mori  <generator> [<generator> ...]
morikrull krull <generator> [<generator> ...]
```

A generator is the minimal polynomial of a positive algebraic number, constant coefficient
first, `p/q` allowed. So `sqrt 2` is `-2,0,1` and the golden ratio is `-1,-1,1`.

```sh
./build/morikrull krull -2,0,1           # sqrt 2
./build/morikrull krull -5,2,-2,1        # the real root of X^3 - 2X^2 + 2X - 5
./build/morikrull krull -5,0,1 -1,-1,1   # sqrt 5 and the golden ratio
```

One generator works in any degree. Several only work when they generate a quadratic field,
which the program finds from the polynomials — you do not have to name it. Append `@` and a
nearby value to pick between two positive roots:

```sh
./build/morikrull krull 5,-5,1@3.6 -5,0,1   # (5 + sqrt 5)/2, not (5 - sqrt 5)/2
```

Output contains the verdict, the reason, and the generating set. For one generator
`alpha^j = [...]` is in the basis `1, alpha, ..., alpha^(d-1)`, and for several, elements print
as `a + b*sqrt(D)` and the set is `G_m` straight from Section 4, so it can carry a monomial
the others already span. A failed Krull test adds the lattice hole witnessing it.

```
$ ./build/morikrull krull -5,2,-2,1
Krull: false
reason: not idv closed in q(H)
additive generators (4):
  alpha^0 = [1, 0, 0]
  alpha^1 = [0, 1, 0]
  alpha^2 = [0, 0, 1]
  alpha^3 = [5, -2, 2]
Hilbert basis (5):
  [0, 0, 1]
  [0, 1, 0]
  [1, 0, 0]
  [3, -1, 1]
  [5, -2, 2]
lattice hole: [3, -1, 1]
```

Bad input exits 1 with a message on stderr.

## Files

    include/morikrull.hpp  library interface
    src/morikrull.cpp      the three procedures of the paper
    src/main.cpp           command line front end
    tests/regression.sh    Examples 4.1, 5.6, 5.7, 5.8, 5.9 and the rejected input
