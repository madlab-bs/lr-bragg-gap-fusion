"""Independent test of the general-N polynomial identities of Section 6.1 (Eq. (6.1)).

Written from the equations, not from the rod or beam code, on random complex hosts:
  "Theorem 1": sin(k1 L) = 0 makes c = C1 a root of P(c) in Eq. (6.1) for every K_T and N;
  "Theorem 2": as K_T -> inf exactly one root diverges, provided sum_j a_j S_j != 0
  (part (b) constructs hosts with sum_j a_j S_j = 0 to show that this condition is needed
  for the branch count; see the discussion after Eq. (6.8)).
The printed labels "Theorem 1/2" refer to these two statements.

Produces the two random-host tests quoted in Section 6.2: the nodal-root identity for N = 1 to 6,
and the count of divergent roots on 2000 random hosts (N = 2 to 6).
Run from python/nchannel/:  python nchan.py   (NumPy; a few seconds; output in nchan_output.txt)
"""
import numpy as np
rng=np.random.default_rng(20260824)
def P(c,C,S,a,KT):
    N=len(C); tot=np.prod(C-c)
    for j in range(N):
        tot+=KT*1j*a[j]*S[j]*np.prod([C[l]-c for l in range(N) if l!=j])
    return tot
def coeffs(C,S,a,KT):
    """polynomial coefficients of P(c) in descending powers"""
    N=len(C); poly=np.poly(C)                       # prod (c - C_j) -> sign (-1)^N * prod(C_j-c)
    poly=poly*((-1)**N)
    for j in range(N):
        others=[C[l] for l in range(N) if l!=j]
        pj=np.poly(others) if others else np.array([1.0])
        pj=pj*((-1)**(N-1))
        poly=poly+KT*1j*a[j]*S[j]*np.concatenate([[0],pj])
    return poly
print("=== Theorem 1: sin(k1 L)=0 => c=C1 is a root, any K_T, any N ===")
worst=0
for N in range(1,7):
    for t in range(300):
        C=rng.normal(size=N)+1j*rng.normal(size=N); S=rng.normal(size=N)+1j*rng.normal(size=N)
        a=rng.normal(size=N)+1j*rng.normal(size=N)
        S[0]=0.0                                     # sin k1 L = 0
        for KT in [1e-3,1.0,1e6,1e12]:
            pc=coeffs(C,S,a,KT); v=abs(np.polyval(pc,C[0]))/max(abs(pc))
            worst=max(worst,v)
print("   max normalized |P(C1)| over N=1..6, 300 random hosts, 4 values of K_T :",f"{worst:.3e}")

print("\n=== Theorem 2: does EXACTLY ONE root diverge as K_T -> inf? ===")
print("   (a) generic hosts")
bad=0; tot=0
for N in range(2,7):
    for t in range(400):
        C=rng.normal(size=N)+1j*rng.normal(size=N); S=rng.normal(size=N)+1j*rng.normal(size=N)
        a=rng.normal(size=N)+1j*rng.normal(size=N)
        r1=np.roots(coeffs(C,S,a,1e8)); r2=np.roots(coeffs(C,S,a,1e14))
        nb1=np.sum(np.abs(r1)>1e4); nb2=np.sum(np.abs(r2)>1e4)
        tot+=1
        if not(nb1==1 and nb2==1): bad+=1
print(f"       {tot-bad}/{tot} cases had exactly one divergent root")
print("   (b) hosts CONSTRUCTED so that sum_j a_j S_j = 0  (the stated hypothesis violated)")
bad=0; tot=0; ex=None
for N in range(2,7):
    for t in range(400):
        C=rng.normal(size=N)+1j*rng.normal(size=N); S=rng.normal(size=N)+1j*rng.normal(size=N)
        a=rng.normal(size=N)+1j*rng.normal(size=N)
        a[-1]=-np.sum(a[:-1]*S[:-1])/S[-1]           # force sum a_j S_j = 0
        r2=np.roots(coeffs(C,S,a,1e14))
        nb=np.sum(np.abs(r2)>1e4); tot+=1
        if nb!=1:
            bad+=1
            if ex is None: ex=(N,nb,np.sort_complex(r2))
print(f"       {tot-bad}/{tot} cases had exactly one divergent root  ->  {bad} VIOLATIONS")
print("       example:",f"N={ex[0]}, divergent roots = {ex[1]} (not 1); roots =",np.array2string(ex[2],precision=3))
