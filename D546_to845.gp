\\ D546_to845.gp -- ICARM curve #845 (rank 31) is a fibre of an elliptic fibration on the member X = surface(pt(2)) of
\\ D546_K3_family.gp.  The fibration is constructed from the family's fibration by two explicit 2-neighbour steps.
\\ usage: gp -q --default parisizemax=2G D546_K3_family.gp D546_to845.gp
\\
\\ Source: https://elliptic-rank.icarm.cloud, curve #845 (curve and points: snapshot of database.json, 2026-09-23).  The
\\ commentary of curve/845.json (by 7fff-zip) calls it "Specialization at T = 504307/1742937 of the rank-16 fibration of
\\ X1092" and gives the model NF845 below, in the family's normal form
\\       y^2 - (L x + B) y = x (x + D)(x + E),   B = p q (L + p + q) - p E - q D.
\\ It is not a member of the 24-I1 family: 96T + 47 divides p and D, and the fibre there is I2 (22 I1 + 1 I2).  So it is a
\\ different elliptic fibration on X (Mordell-Weil rank 16, plus the A1 of the I2 fibre).  X is also the surface of #302
\\ (D546_to302.gp).
\\
\\ Notation.  X has fibre class F, zero section O, sections P_1..P_17 = surface(pt(2))[3], and 24 fibres I1.
\\ Step 0.  The construction starts from X after the change of coordinates v0 below (a Moebius change of T and a scaling).
\\ NS(X) = U + W(-1), where W is the Mordell-Weil lattice; <P,Q> = 2 + P.O + Q.O - P.Q, and h(P) = 4 + 2 P.O.
\\
\\ Step 1 (rank 17 -> rank 17).  S = sum W1_i P_i has h(S) = 12, i.e. S.O = 4.  The class F' = O + S - F has
\\   F'^2 = 0, F'.F = 2, F'.O = F'.S = 1: an elliptic fibration with sections O and S.  W1 is minimal in W1 + 2W,
\\   so F' is nef and the new fibres are again irreducible.  Elliptic parameter u = u1/u2, (u1, u2) a basis of
\\   H^0(O + S - F):  u = (A0 + B0 Z lam)/(A1 + B1 Z lam),  lam = (y + yS)/(x - xS),  xS = XS/Z^2 (deg Z = 4).
\\   Fixing u turns the curve into w^2 = quartic in T over Q(u); O gives a rational point, and Connell's formulas a
\\   Weierstrass model.  An old section P with P.F' = 1 has u|P of degree 1 in T and becomes a section tau(P).
\\   The model is then moved by u = (-11T - 12)/(-4T + 3) and scaled by 3^8.
\\ Step 2 (rank 17 -> rank 16).  S2 is a section of step 1 with S2.O = 2.  F'' = O + S2 has F''^2 = 0, F''.F' = 2,
\\   and the curve O + S2 is itself a fibre: an I2 at u = oo.  O is a fibre component, so the zero section is a
\\   section Pb with Pb.F'' = 1 (Pb.O = 0, <Pb, S2> = 3).  H^0(O + S2) = <1, u>, u = (A0 + B0 Z lam)/Z^2 (deg Z = 2).
\\   Result: 22 I1 + I2, Mordell-Weil rank 16.  S2 and Pb are written below as sums of sections tau(P_v) of step 1.
\\   The Weierstrass model is the Jacobian of the fibration and does not depend on which section is the zero
\\   section, so any Pb with Pb.F'' = 1 gives the same model.
\\ Step 3.  u = (22 - 84T)/(96T + 47) (the I2 goes to T = -47/96) and a rational scaling give the #845 model NF845
\\   over Q(T).  Its fibre at T0 = 504307/1742937 is #845.
\\ The classes W1, S2TERMS, VB (in the basis P_1..P_17) and the maps v0, M1, M2 were found by a search; this file only
\\ checks them.
\\
\\ What this file checks: the three steps above, the fibres of the #845 model, that its fibre at T0 is #845, and that
\\ the 31 published points lie on #845 and are independent.  Sections of the rank-16 fibration over Q(T) are not given.

N_U = varlower("u");
\\ ---- curve y^2 = x^3 + C[1] x^2 + C[2] x + C[3] over Q(T); points [x, y], [0] = O ----
n_cs(ai) = [ai[2] + ai[1]^2/4, ai[4] + ai[1]*ai[3]/2, ai[5] + ai[3]^2/4];
n_cspt(ai, P) = [P[1], P[2] + (ai[1]*P[1] + ai[3])/2];
n_add(C, P, Q) = {
  my(l, x3);
  if(P == [0], return(Q)); if(Q == [0], return(P));
  if(P[1] == Q[1], if(P[2] == -Q[2], return([0])); l = (3*P[1]^2 + 2*C[1]*P[1] + C[2])/(2*P[2]), l = (Q[2] - P[2])/(Q[1] - P[1]));
  x3 = l^2 - C[1] - P[1] - Q[1]; [x3, -(P[2] + l*(x3 - P[1]))];
}
n_neg(P) = if(P == [0], P, [P[1], -P[2]]);
n_mul(C, P, n) = { my(R = [0], Q = if(n < 0, n_neg(P), P)); for(i = 1, abs(n), R = n_add(C, R, Q)); R; }
n_lin(C, pts, c) = { my(R = [0]); for(i = 1, #c, if(c[i], R = n_add(C, R, n_mul(C, pts[i], c[i])))); R; }
n_on(C, P) = P[2]^2 == P[1]^3 + C[1]*P[1]^2 + C[2]*P[1] + C[3];
n_dotO(P) = { my(n = numerator(P[1]), d = denominator(P[1])); poldegree(d, 'T)/2 + max(0, poldegree(n, 'T) - poldegree(d, 'T) - 4)/2; }
n_ht(P) = 4 + 2*n_dotO(P);
n_pair(C, P, Q) = { my(R = n_add(C, P, n_neg(Q))); (n_ht(P) + n_ht(Q) - if(R == [0], 0, n_ht(R)))/2; }
\\ c4, c6 in Q(u) -> [c4n, c6n, r]: c4 = r^4 c4n, c6 = r^6 c6n, c4n, c6n in Z[u] with no p^4 | c4n, p^6 | c6n
n_norm(c4, c6) = {
  my(ps = [], r = 1, c4n, c6n, den, k, g4, g6, g);
  foreach([c4, c6], c, if(c, foreach([numerator(c), denominator(c)], f, if(poldegree(f, N_U) > 0, ps = setunion(ps, Set(factor(f)[,1]~))))));
  foreach(ps, pi, r *= pi^min(floor(if(c4, valuation(c4, pi), 10^9)/4), floor(if(c6, valuation(c6, pi), 10^9)/6)));
  c4n = c4/r^4; c6n = c6/r^6; c4n = numerator(c4n)/denominator(c4n); c6n = numerator(c6n)/denominator(c6n);
  den = lcm(concat([denominator(z) | z <- Vec(c4n)], [denominator(z) | z <- Vec(c6n)])); k = 1/den; c4n /= k^4; c6n /= k^6;
  g4 = content(c4n); g6 = content(c6n); g = gcd(g4, g6);
  if(g != 1, my(f = factor(g)); for(i = 1, #f~, my(p = f[i,1], m = min(valuation(g4, p)\4, valuation(g6, p)\6));
    if(m > 0, c4n /= p^(4*m); c6n /= p^(6*m); k *= p^m)));
  [c4n, c6n, r*k];
}
\\ Connell: w^2 = e[5] s^4 + e[4] s^3 + e[3] s^2 + e[2] s + alpha^2 (alpha^2 = e[1]) -> [a1, a2, a3, a4, a6]
n_connell(e, alpha) = { my(a2 = e[3] - e[2]^2/(4*alpha^2), a4 = -4*alpha^2*e[5]); [e[2]/alpha, a2, 2*alpha*e[4], a4, a2*a4]; }
n_c4c6(a) = { my(b2 = a[1]^2 + 4*a[2], b4 = a[1]*a[3] + 2*a[4], b6 = a[3]^2 + 4*a[5]); [b2^2 - 24*b4, -b2^3 + 36*b2*b4 - 216*b6, b2]; }
\\ the point (s, w) of the quartic -> the point of Y^2 = X^3 - 27 c4n X - 54 c6n
n_qpt(qc, s, w) = {
  my([alpha, e, a, b2, r] = qc, X, Y);
  if(s == 0, X = -a[2]; Y = a[1]*a[2] - a[3],
    X = (2*alpha*(w + alpha) + e[2]*s)/s^2;
    Y = (4*alpha^2*(w + alpha) + 2*alpha*(e[2]*s + e[3]*s^2) - e[2]^2*s^2/(2*alpha))/s^3);
  [(36*X + 3*b2)/r^2, 108*(2*Y + a[1]*X + a[3])/r^3];
}
\\ the fibre over u: lam = N/(Z l), and w = 2x + xS + A2 - lam^2 satisfies w^2 = Z^2 q(T)/l^4, q = P/Z^6
n_quartic(C, XS, YS, Z, N, l) = {
  my(P = N^4 - 2*(3*XS + C[1]*Z^2)*N^2*l^2 - 8*YS*N*l^3 + (C[1]^2*Z^4 - 2*C[1]*XS*Z^2 - 3*XS^2 - 4*C[2]*Z^4)*l^4, q = P/Z^6);
  if(denominator(q) != 1 || poldegree(q, 'T) != 4, error("the quartic has the wrong shape")); q;
}
n_root1(f) = -polcoef(f, 0, 'T)/polcoef(f, 1, 'T);       \\ the root in Q(u) of f, of degree 1 in T
n_Qmodel(q, t0, alpha) = {                                \\ Connell data of the quartic at its point over T = t0
  my(e = vector(5, i, polcoef(subst(q, 'T, t0 + 'T), i-1, 'T)), a, cc, res);
  if(e[1] != alpha^2, error("the base point is not on the quartic"));
  a = n_connell(e, alpha); cc = n_c4c6(a); res = n_norm(cc[1], cc[2]);
  [res[1], res[2], [alpha, e, a, cc[3], res[3]]];
}

\\ ---- step 1: 2-neighbour F' = O + S - F (S.O = 4).  Returns [c4n, c6n, ctx], c4n, c6n in Z[u] ----
nb17(C, S) = {
  my(Z = denominator(S[1]), XS, YS, cols, K, A0, B0, A1, B1, l, N, tO, R);
  if(!issquare(Z, &Z) || poldegree(Z) != 4, error("S.O != 4"));
  XS = S[1]*Z^2; YS = S[2]*Z^3;
  \\ H^0(O + S - F): (A0 + B0 Z lam)/Z^2 with deg A0 <= 7, deg B0 <= 1, A0 XS = B0 YS mod Z^2
  cols = concat(vector(8, k, ('T^(k-1)*XS) % Z^2), vector(2, j, (-'T^(j-1)*YS) % Z^2));
  K = matkerint(matrix(8, 10, i, j, polcoef(cols[j], i-1, 'T))); if(#K != 2, error("h^0 != 2"));
  A0 = Pol(Vecrev(K[1..8,1]), 'T); B0 = Pol(Vecrev(K[9..10,1]), 'T); A1 = Pol(Vecrev(K[1..8,2]), 'T); B1 = Pol(Vecrev(K[9..10,2]), 'T);
  N = A0 - N_U*A1; l = N_U*B1 - B0;
  tO = n_root1(l);                                         \\ O and S lie over T = tO
  R = n_Qmodel(n_quartic(C, XS, YS, Z, N, l), tO, subst(N, 'T, tO)^2/subst(Z, 'T, tO)^3);
  [R[1], R[2], [C, S, Z, A0, B0, A1, B1, tO, R[3]]];
}
tau17(ctx, P) = {                                          \\ the old section P (P.F' = 1) on the step-1 model
  my([C, S, Z, A0, B0, A1, B1, tO, qc] = ctx, lam, uP, t, lt);
  if(P == S, return(n_qpt(qc, 0, 0)));                     \\ S: the second point over T = tO
  lam = (P[2] + S[2])/(P[1] - S[1]);
  uP = (A0 + B0*Z*lam)/(A1 + B1*Z*lam);
  if(max(poldegree(numerator(uP), 'T), poldegree(denominator(uP), 'T)) != 1, error("P.F' != 1"));
  t = n_root1(numerator(uP) - N_U*denominator(uP)); lt = subst(lam, 'T, t);
  n_qpt(qc, t - tO, (2*subst(P[1], 'T, t) + subst(S[1] + C[1], 'T, t) - lt^2)*(N_U*subst(B1, 'T, t) - subst(B0, 'T, t))^2/subst(Z, 'T, t));
}
\\ ---- step 2: 2-neighbour F'' = O + S (S.O = 2), zero section Pb.  Returns [c4n, c6n] ----
nb16(C, S, Pb) = {
  my(Z = denominator(S[1]), XS, YS, cols, K, j, A0, B0, lam, uP, t, lt);
  if(!issquare(Z, &Z) || poldegree(Z) != 2, error("S.O != 2"));
  XS = S[1]*Z^2; YS = S[2]*Z^3;
  \\ H^0(O + S) = <1, u>, u = (A0 + B0 Z lam)/Z^2, deg A0 <= 4, B0 constant
  cols = concat(vector(5, k, ('T^(k-1)*XS) % Z^2), [(-YS) % Z^2]);
  K = matkerint(matrix(4, 6, i, j, polcoef(cols[j], i-1, 'T))); if(#K != 2, error("h^0 != 2"));
  j = if(K[6,1], 1, 2); A0 = Pol(Vecrev(K[1..5,j]), 'T); B0 = K[6,j];
  lam = (Pb[2] + S[2])/(Pb[1] - S[1]); uP = (A0 + B0*Z*lam)/Z^2;
  if(max(poldegree(numerator(uP), 'T), poldegree(denominator(uP), 'T)) != 1, error("Pb.F'' != 1"));
  t = n_root1(numerator(uP) - N_U*denominator(uP)); lt = subst(lam, 'T, t);   \\ Pb meets the fibre over u at T = t
  n_Qmodel(n_quartic(C, XS, YS, Z, N_U*Z^2 - A0, B0), t,
           (2*subst(Pb[1], 'T, t) + subst(S[1] + C[1], 'T, t) - lt^2)*B0^2/subst(Z, 'T, t))[1..2];
}
\\ f(u) of weight k -> f((a T + b)/(c T + d)) (c T + d)^k
n_mob(f, m, k) = subst(f, N_U, (m[1,1]*'T + m[1,2])/(m[2,1]*'T + m[2,2]))*(m[2,1]*'T + m[2,2])^k;

\\ ============================================= data =============================================
\\ step 0: the frame of the construction.  f of weight k (1: L,p,q;  2: D,E,x;  3: y) -> s^k (32T + 5)^(2k) f((76T - 5)/(32T + 5)),
\\ s = 3^4/2^3: a Moebius change of T and a scaling, applied to surface(pt(2)) and its 17 sections.
v0(f, k) = (3^4/2^3)^k * (32*'T + 5)^(2*k) * subst(f, 'T, (76*'T - 5)/(32*'T + 5));
W1 = [0, 2, 0, 2, -1, -1, 0, 1, -2, 1, -2, -2, 3, -1, -1, 1, 0];       \\ step 1: S = sum W1_i P_i, h(S) = 12
M1 = [-11, -12; -4, 3]; MU1 = 3^8;                                      \\ step 1: Moebius map and scaling
\\ step 2: S2 = S2_CS tau(S) + sum c tau(P_v) over S2TERMS = [c, v];  Pb = -tau(P_VB)
S2_CS = -1;
{S2TERMS = [
  [-1, [1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1]],
  [ 1, [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0]],
  [ 1, [0, 0, 0, 1, -1, 0, 0, 0, -1, 1, 0, -1, 1, 0, 0, 0, 0]],
  [ 1, [-1, 0, 0, 1, 0, 0, 0, 0, -1, 1, -1, -1, 1, -1, -1, 0, 0]],
  [-1, [0, 0, 0, 1, 0, 0, 0, 0, -1, 1, -1, -1, 1, 0, -1, 0, 0]],
  [-1, [-1, 0, -1, 1, 0, 0, 0, 0, -1, 1, -1, -1, 1, -1, -1, 0, 0]],
  [ 1, [0, 0, -1, 1, 1, 0, 0, 0, 0, 1, -1, -1, 1, 0, -1, 0, 1]],
  [ 1, [0, 1, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 1]]
];}
VB = [1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1];
M2 = [-84, 22; 96, 47];                                                 \\ step 3: u = (22 - 84 T)/(96 T + 47)
T0 = 504307/1742937;
\\ #845 (ICARM) and the model NF845 of its commentary, [L, p, q, D, E]
E845 = [1, 0, 0, -68906458325887563494253912743860596601427440743924636890868535, 219726256542658179767501421714589587209898779390156657938983977041099527263435646417385213097];
{NF845 = [102579386*'T^2 - 110801088*'T - 18001125, -(96*'T + 47)*(1094441*'T - 760632), 101*(163604*'T^2 + 454821*'T - 85571),
  -(96*'T + 47)*(7477016504034*'T^3 + 7131320229659*'T^2 - 30633770734995*'T + 6299022506152),
  -101*(58579461360*'T^4 + 18954357876307*'T^3 - 8762256702014*'T^2 - 15609000269917*'T + 2960508883604)];}

{PTS845 = [
  [4492189723065914919256867882934, 914724714768630901238259806376949142963232533],
  [4567121567039911266825199317014, 534643646961232792435902595697216090643639413],
  [4586421283438238075388908357654, 410821791963280157575133208491598963338280053],
  [5037509468610356621679378871334, 666098792397671427124894554496634456754265733],
  [2659630394938362489510492764438, 7434634969447538378141678087206140041567980661],
  [3663881909579733757226379828044, 4055267373636724119245464691819620518768714863],
  [3677812647428541656906378144054, 4006049353458969139310786882901444065827897493],
  [4098341294835689182470629375254, 2482228482598568995384324738277746299302465653],
  [4164965728556382418786791819254, 2232180041009927043930804781152182748288551093],
  [4210076388396811999333015844534, 2060904208652281762820569092863917211521306133],
  [4417431261416015121795299132054, 1239813367034742170365682325542949695429486453],
  [4431929075218277083918138761814, 1178834426993205829246391686760443853789620213],
  [5105676126910662784165003206038, 1003293868066583439325644633458295443785691381],
  [5176341638430813588645978697934, 1319257072366563187890692056699670017335927533],
  [5316559979035457788289474404022, 1912540767664534490942467385279679856808556053],
  [5324224606886431739743202668022, 1944363137574878874350541396985641121979088053],
  [5482287408296083204197995043734, 2594995395099938359808519242326412410251570933],
  [5548697434495430663915305345814, 2866828288411478685823776455429234788876484213],
  [5615878802954628697457810529302, 3141664450647791215948876141974015871913938805],
  [5733477463417257660942569682134, 3623241806674604699060111526578454755278054133],
  [5874686722234872651737487188438, 4203481938879821615647418800850411741692040501],
  [7017936048521021994702580091414, 9043699426822951025499785852469862997523101813],
  [757752158220185258826363604934, 12959449505271470233127111406483601080657066533],
  [973557344004335969787009680534, 12392119236646578470308138231491128734565177333],
  [1013068729907318835646725610934, 12286537312070329678298138251534690134629856533],
  [1639444315586958180395107629614, 10543453748450346001731175965059489249608951373],
  [226164351008217133335210359392886/49, 36401096114268784547188426613387846966845341059/343],
  [123713595815943158865232342397366/49, 2697460873874820084336950933813275583929181950019/343],
  [3708780681953243884687375900371494/841, 30996726003060644535638723469585545194108289295337/24389],
  [-1732275322545229786628751579745114/1369, 884459401907288495674534376294947971656404061298689/50653],
  [3228375666134932050157146796346896575014/650199001, 159093427395242901986547740315901913477350720665857393087/16579424326499]];}

\\ ============================================= run =============================================
n_chk(s, c) = { printf("  %-100s %s\n", s, if(c, "ok", "FAILED")); c; }
{
  my(t0 = getabstime(), ok = 1, H, X, A, C0, P0, S, R1, c4, c6, C8, tr, S2, Pb, R2, E, r4, r6, lam2);
  print("D546_to845: #845 from surface(pt(2)) of D546_K3_family.gp");
  X = surface(pt(2)); NF0 = vector(5, i, v0(X[1][i], if(i < 4, 1, 2))); A = ainvs(NF0);
  P0 = [[v0(P[1], 2), v0(P[2], 3)] | P <- X[3]];
  ok = n_chk("step 0: surface(pt(2)) and its 17 sections in the frame of the construction", vecmin(vector(#P0, i, oncurve(A, P0[i])))) && ok;
  C0 = n_cs(A); P0 = [n_cspt(A, P) | P <- P0];
  S = n_lin(C0, P0, W1);
  ok = n_chk(Str("step 1: S = sum W1_i P_i has S.O = ", n_dotO(S), " (height ", n_ht(S), ")"), n_dotO(S) == 4) && ok;
  R1 = nb17(C0, S);
  c4 = n_mob(R1[1], M1, 8)/MU1^4; c6 = n_mob(R1[2], M1, 12)/MU1^6; C8 = [0, -27*c4, -54*c6];
  ok = n_chk(Str("step 1: Y^2 = X^3 - 27 c4 X - 54 c6, c4, c6 in Z[T] of degrees ", [poldegree(c4), poldegree(c6)], ", Delta of degree 24, squarefree"),
             denominator(content([c4, c6])) == 1 && poldegree(c4^3 - c6^2) == 24 && poldisc(c4^3 - c6^2) != 0) && ok;
  tr = (P -> my(Q = tau17(R1[3], P)); [n_mob(Q[1], M1, 4)/MU1^2, n_mob(Q[2], M1, 6)/MU1^3]);
  S2 = n_mul(C8, tr(S), S2_CS);
  foreach(S2TERMS, z, S2 = n_add(C8, S2, n_mul(C8, tr(n_lin(C0, P0, z[2])), z[1])));
  Pb = n_neg(tr(n_lin(C0, P0, VB)));
  ok = n_chk("step 2: S2 and Pb lie on the step-1 model", n_on(C8, S2) && n_on(C8, Pb)) && ok;
  ok = n_chk(Str("step 2: S2.O = ", n_dotO(S2), ", Pb.O = ", n_dotO(Pb), ", <Pb, S2> = ", n_pair(C8, Pb, S2), "   (F'' = O + S2, Pb.F'' = 1)"),
             n_dotO(S2) == 2 && n_dotO(Pb) == 0 && n_pair(C8, Pb, S2) == 3) && ok;
  R2 = nb16(C8, S2, Pb);
  ok = n_chk(Str("step 2: rank-16 model over u: deg c4 = ", poldegree(R2[1]), ", Delta of degree ", poldegree(R2[1]^3 - R2[2]^2),
                 ", squarefree: 22 I1 + I2 at u = oo"), poldegree(R2[1]^3 - R2[2]^2) == 22 && poldisc(R2[1]^3 - R2[2]^2) != 0) && ok;
  E = ellinit(ainvs(NF845));
  r4 = n_mob(R2[1], M2, 8)/E.c4; r6 = n_mob(R2[2], M2, 12)/E.c6;
  lam2 = if(poldegree(r4, 'T) == 0 && poldegree(r6, 'T) == 0 && r4, r6/r4, 0);
  ok = n_chk(Str("step 3: after u = (22 - 84T)/(96T + 47): c4 = l^4 c4(#845 model), c6 = l^6 c6(#845 model), l^2 = ", lam2),
             lam2 && r4 == lam2^2 && r6 == lam2^3) && ok;
  ok = n_chk("step 3: l is rational, so the two models are isomorphic over Q(T) (no quadratic twist)", issquare(lam2)) && ok;
  ok = n_chk("the fibre of the #845 model at T0 = 504307/1742937 has the minimal model of #845",
             ellminimalmodel(ellinit(subst(ainvs(NF845), 'T, T0)))[1..5] == ellminimalmodel(ellinit(E845))[1..5]) && ok;
  ok = n_chk(Str("fibres of the #845 model: ", fibres(ainvs(NF845)), ", the I2 at 96T + 47 = 0"),
             fibres(ainvs(NF845)) == [["I1", 22], ["I2", 1]] && poldegree(gcd(disc_(ainvs(NF845)), (96*'T + 47)^2), 'T) == 2) && ok;
  H = ellheightmatrix(ellinit(E845), PTS845);
  ok = n_chk(Str("the ", #PTS845, " published points lie on #845; rank of their height matrix: ", matrank(H)),
             vecmin(vector(#PTS845, i, ellisoncurve(ellinit(E845), PTS845[i]))) && matrank(H) == #PTS845) && ok;
  printf("D546_to845: %s   (%.1f s)\n", if(ok, "all checks passed", "SOME CHECK FAILED"), (getabstime() - t0)/1000.);
}
