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

Without Nix, put any recent `gp` on your PATH and run `./check.sh`. All checks together take about 15 seconds.

For each family, `check.gp` runs:

1. `verify_generic()`: exact computations over the function field Q(C). These give statements (1)–(4) of the
   theorem: every section lies on the surface; Δ has degree 24 and is squarefree; the Gram matrix has the stated
   determinant and discriminant group, with no roots and no even overlattice; and the family is not isotrivial.
2. `verify(P)` at every point P = mG + T with |m| ≤ 8 and T ∈ C[2] (68 points). This certifies each member directly.
   The points that fail are exactly the degenerate ones listed in the header: m ∈ {−1, 0, 1} for D510 and
   m ∈ {−2, −1, 0, 1} for D546.

`check.sh` also runs `D546_to302.gp` and `D546_to845.gp` (see below). It writes the output to `out/`, removes the timings, and compares the result with `expected/`.

For interactive use:

    gp -q D546_K3_family.gp        # start a fresh session: the file fixes the variable order T > Y > X
    ? M = surface(pt(6));          # 6G -> [[L,p,q,D,E,B], [a1,...,a6], 17 sections]
    ? verify(pt(6))

## Relation to ICARM curves #302 and #845

Both ICARM curves are fibres of elliptic fibrations on the same surface: the D546 member at P = 2G. The curve data are
from a snapshot of the [ICARM elliptic-curve rank database](https://elliptic-rank.icarm.cloud) taken on 2026-09-23.
Each file is read after the family file, e.g. `gp -q D546_K3_family.gp D546_to302.gp`.

- `D546_to302.gp`: #302 (rank 31) is a fibre of the family's own fibration. A change of T and a scaling map
  `surface(pt(2))` exactly onto the model published in the #302 commentary; #302 is its fibre at
  T = 164518/924945. The 17 sections carry over, so that model has Mordell–Weil rank 17 over Q(T).
- `D546_to845.gp`: #845 (rank 31) is a fibre of a different fibration, with 22 I₁ + I₂ fibres and Mordell–Weil
  rank 16. The file constructs it from `surface(pt(2))` by two explicit 2-neighbour steps. The result is isomorphic
  over Q(T) to the model in the #845 commentary, and #845 is its fibre at T = 504307/1742937. The sections of the
  rank-16 fibration are not given.

Both files also check that the 31 published points lie on the curve and are independent.

## Input from the literature

The proof uses two results that it quotes rather than computes:

- H. Sterk, *Finiteness results for algebraic K3 surfaces* (1985). A K3 surface has only finitely many elliptic
  fibrations up to automorphisms, so a non-isotrivial family has infinitely many non-isomorphic members.
- M. Schütt, *Fields of definition of singular K3 surfaces* (2007). This bounds the rank over Q when X_P has
  geometric Picard number 20.
