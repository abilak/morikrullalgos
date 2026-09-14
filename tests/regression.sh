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

expect 'Mori: true' quadratic-mori 2 0,1
expect 'Mori: false' quadratic-mori 2 2,1
expect 'Krull: true' quadratic-krull 2 0,1
expect 'Krull: false' quadratic-krull 2 0,2 0,3
expect 'lattice hole: 0 + 1*sqrt(2)' quadratic-krull 2 0,2 0,3
expect 'Krull: true' quadratic-krull 5 0,1 1/2,1/2

expect 'Mori: true' cyclic-mori -2,0,1
expect 'Mori: false' cyclic-mori 1,-3,1
expect 'Mori: false' cyclic-mori -1/2,1
expect 'Krull: true' cyclic-krull -2,0,1
expect 'Krull: true' cyclic-krull -2,0,0,1
expect 'Krull: false' cyclic-krull -5,2,-2,1
expect 'lattice hole: [3, -1, 1]' cyclic-krull -5,2,-2,1

printf 'all regression cases passed\n'
