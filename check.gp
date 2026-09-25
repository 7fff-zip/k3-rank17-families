\\ check.gp -- read after a family file:  gp -q D510_K3_family.gp check.gp
\\ 1. verify_generic(): the exact computations over the function field of C behind the THEOREM in the file header.
\\ 2. verify(P) at every point P = mG + T, |m| <= 8, T in C[2]; lists the points that are not certified.

ok = verify_generic();
print();
bad = List();
{for(m = -8, 8, for(t = 0, 3,
  r = iferr(verify(pt(m, t), []), E, printf("mG + T for m = %d, t = %d: %s\n", m, t, errname(E)); 0);
  if(!r, listput(bad, [m, t]))))}
printf("\nnot certified [m, t]: %s\n", Vec(bad));
printf("check: %s\n", if(ok, "verify_generic passed", "verify_generic FAILED"));
\q
