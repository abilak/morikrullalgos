#!/bin/sh
set -eu

program=$1

expect() {
  wanted=$1
  shift
  output=$($program "$@")
  printf '%s\n' "$output" | grep -F "$wanted" >/dev/null || {
    printf 'expected %s in output of: %s\n%s\n' "$wanted" "$*" "$output" >&2
    exit 1
  }
}

expect_error() {
  wanted=$1
  shift
  if output=$($program "$@" 2>&1); then
    printf 'expected a failure from: %s\n%s\n' "$*" "$output" >&2
    exit 1
  fi
  printf '%s\n' "$output" | grep -F "$wanted" >/dev/null || {
    printf 'expected error %s from: %s\n%s\n' "$wanted" "$*" "$output" >&2
    exit 1
  }
}

# cor 3.13, one generator in any degree
expect 'Mori: true'  mori -2,0,1
expect 'Mori: true'  mori -1,-1,1
expect 'Mori: false' mori 1,-3,1
expect 'Mori: false' mori -1,1,0,1
expect 'Mori: false' mori -1/2,1
expect 'Krull: true' krull -2,0,1
expect 'Krull: true' krull -1,-1,1
expect 'Krull: true' krull -2,0,0,1

# example 5.6, d = 3 and k = 4
expect 'Krull: false' krull -5,2,-2,1
expect 'alpha^3 = [5, -2, 2]' krull -5,2,-2,1
expect 'lattice hole: [3, -1, 1]' krull -5,2,-2,1

# thm 4.2, several generators sharing a quadratic field
expect 'Mori: true' mori -2,0,1 -3,1
expect 'Mori: true' mori -5,0,1 -1,-1,1
expect 'Krull: true' krull -5,0,1 -1,-1,1                 # example 5.9
expect 'is not an algebraic integer' mori -2,0,1 -7,-4,4  # example 4.1
expect 'no generator has a negative conjugate' mori 2,-4,1@3.4 1,-6,1@5.8

# example 5.8, 2sqrt2 and 3sqrt2
expect 'Mori: true' mori -8,0,1 -18,0,1
expect 'Krull: false' krull -8,0,1 -18,0,1
expect 'lattice hole: 0 + 1*sqrt(2)' krull -8,0,1 -18,0,1

# the field comes out of the polynomials, and @ picks between two positive roots
expect '5/2 + 1/2*sqrt(5)' krull 5,-5,1@3.6 -5,0,1
expect 'is not weak Perron' krull 5,-5,1@1.4 -5,0,1

# rejected input
expect_error 'reducible' mori 6,-5,1
expect_error 'reducible' mori 0,0,1
expect_error 'degree at least one' mori 3
expect_error '@ does nothing' mori -2,0,1@1.4
expect_error 'different quadratic fields' krull -2,0,1 -5,0,1
expect_error 'span Q' krull -3,1 -5,1
expect_error 'degree 1 or 2' krull -2,0,1 -2,0,0,1
expect_error 'rational roots' krull -2,0,1 4,-4,1
expect_error 'not real' krull -2,0,1 1,0,1
expect_error 'must be positive' krull -2,0,1 -5,-4,4@-0.6
expect_error 'two positive roots' krull 5,-5,1 -5,0,1

printf 'all tests passed\n'
