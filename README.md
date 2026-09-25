# Two one-parameter families of rank-17 elliptic K3 surfaces

This repository contains explicit families of elliptic K3 surfaces over Q(T) with Mordell–Weil rank 17. Each family
is parametrised by the rational points of a rank-1 elliptic curve C, so it has infinitely many members that are
pairwise non-isomorphic over Q̄.

| file                | parameter curve C                          | conductor | generator G | MW lattice                    |
|---------------------|--------------------------------------------|-----------|-------------|-------------------------------|
| `D510_K3_family.gp` | y² + xy + y = x³ + x² − 421x − 3157        | 510       | (−13, 24)   | det 1020, group Z/1020, no roots |
| `D546_K3_family.gp` | y² + xy + y = x³ − 137x + 380              | 546       | (−6, 34)    | det 1092, group Z/1092, no roots |

For P = (X, Y) on C, the member X_P is

    y² − (L x + B) y = x (x + D)(x + E),     B = pq(L + p + q) − pE − qD,

where L, p and q have degree 2 in T, and D and E have degree 4. Their coefficients lie in Z[X, Y]. The 17 sections are
(0,0), (−D,0), (−pq, q(pq − D)) and 14 sections on lines through them. Each `.gp` file is self-contained. Its header
states the theorem, gives the proof and lists usage; the heuristic facts it records are marked as not proven.

## Reproducing

You need PARI/GP (tested with 2.17.3). `flake.nix` and `flake.lock` pin this version:

    nix develop -c ./check.sh

Without Nix, put any recent `gp` on your PATH and run `./check.sh`. Both families together take about 15 seconds.

For each family, `check.gp` runs:

1. `verify_generic()`: exact computations over the function field Q(C). These give statements (1)–(4) of the
   theorem: every section lies on the surface; Δ has degree 24 and is squarefree; the Gram matrix has the stated
   determinant and discriminant group, with no roots and no even overlattice; and the family is not isotrivial.
2. `verify(P)` at every point P = mG + T with |m| ≤ 8 and T ∈ C[2] (68 points). This certifies each member directly.
   The points that fail are exactly the degenerate ones listed in the header: m ∈ {−1, 0, 1} for D510 and
   m ∈ {−2, −1, 0, 1} for D546.

`check.sh` writes the output to `out/`, removes the timings, and compares the result with `expected/`.

For interactive use:

    gp -q D546_K3_family.gp        # start a fresh session: the file fixes the variable order T > Y > X
    ? M = surface(pt(6));          # 6G -> [[L,p,q,D,E,B], [a1,...,a6], 17 sections]
    ? verify(pt(6))

## Input from the literature

The proof uses two results that it quotes rather than computes:

- H. Sterk, *Finiteness results for algebraic K3 surfaces* (1985). A K3 surface has only finitely many elliptic
  fibrations up to automorphisms, so a non-isotrivial family has infinitely many non-isomorphic members.
- M. Schütt, *Fields of definition of singular K3 surfaces* (2007). This bounds the rank over Q when X_P has
  geometric Picard number 20.
