"""2D sweep of the resonator tuning w_r through coincidence with the Bragg frequency w_B.

Exploratory lattice-Green-function calculation for square resonator lattices on a membrane
and a Kirchhoff plate (manuscript Section 7.3).  It is not part of the one-dimensional proof
and is not used as evidence for the criterion.

Produces the sweep summarized in Section 7.3: membrane and plate, contact radius 0.1 L,
mass ratios 1 and 4, tunings 0.90-1.10 w_B, 18 reciprocal-lattice vectors in each direction,
600 Bloch phases, band edges bisected to 1e-8 in w.
Run from python/exploratory_2d/:  python coincidence_sweep.py
(NumPy and SciPy; output in coincidence_sweep_output.txt; the per-test timings vary by machine)

Two design choices:

(1) Robust edge classification -- the omega grid is ANCHORED on w_B exactly, and the "edge
versus interior" question is answered by testing propagation at w_B +- delta over several
decades of delta, never by comparing a bisected edge against w_B with a fixed tolerance.

(2) Speed:
  (a) pole crossings: |q0+G| is smooth in mu, so per reciprocal vector only its min and max
      over the path can bracket k.  Precompute those once; for a given k only the handful of
      columns with min <= k <= max can cross, instead of an argmin over all ~3000 columns.
  (b) the structural sum is evaluated on a coarse mu grid during the omega scan and a fine one
      only when refining an edge.
"""
import numpy as np
from scipy.special import j1

L = 1.0; A = L*L

class Path:
    def __init__(self, kind, nmu=600, acont=0.10, nmax=18):
        self.kind = kind
        self.mux = np.clip(np.linspace(0, np.pi, nmu), 1e-7, np.pi - 1e-10)
        p = np.arange(-nmax, nmax+1)
        P, Q = np.meshgrid(p, p, indexing='ij')
        Gx = (2*np.pi/L)*P.ravel(); Gy = (2*np.pi/L)*Q.ravel()
        qx = self.mux[:, None]/L + Gx[None, :]; qy = 0.0 + Gy[None, :]
        self.qm = np.sqrt(qx*qx + qy*qy); self.q2 = self.qm**2
        x = self.qm*acont*L
        self.F2 = np.where(x > 1e-12,
                           (2*j1(np.where(x > 1e-12, x, 1.0))/np.where(x > 1e-12, x, 1.0))**2, 1.0)
        self.cmin = self.qm.min(axis=0); self.cmax = self.qm.max(axis=0)   # speed-up (2)(a)
    def k(self, w):  return w if self.kind == 'membrane' else np.sqrt(w)
    def S(self, w):
        k = self.k(w)
        return (np.sum(self.F2/(self.q2 - k*k), axis=1)/A if self.kind == 'membrane'
                else np.sum(self.F2/(self.q2*self.q2 - k**4), axis=1)/A)
    def cuts(self, w):
        k = self.k(w)
        sel = np.where((self.cmin <= k) & (self.cmax >= k))[0]
        if sel.size == 0: return np.array([], dtype=int)
        d = self.qm[:, sel] - k
        ch = np.where(np.any(np.sign(d[:-1, :]) != np.sign(d[1:, :]), axis=1))[0]
        return ch

def KT(w, wr, mr):
    K = mr*wr*wr
    return K*mr*w*w/(mr*w*w - K)

def prop(path, w, wr, mr):
    if abs(w - wr) < 1e-13*max(wr, 1.0): return False
    Phi = 1.0 + KT(w, wr, mr)*path.S(w)
    c = np.concatenate(([0], path.cuts(w)+1, [len(Phi)]))
    for i in range(len(c)-1):
        seg = Phi[c[i]:c[i+1]]; seg = seg[np.isfinite(seg)]
        if seg.size >= 2 and seg.min() < 0 < seg.max(): return True
    return False

def edge_or_interior(path, wB, wr, mr, decades=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)):
    """Classify w_B by testing propagation on both sides over several decades (choice (1))."""
    below = [prop(path, wB*(1-d), wr, mr) for d in decades]
    above = [prop(path, wB*(1+d), wr, mr) for d in decades]
    return below, above

def passband_around(path, wB, wr, mr, half=0.45, n=500, tol=1e-8):
    """the propagating interval containing or adjacent to w_B, found on a grid anchored on w_B"""
    lo_w, hi_w = wB*(1-half), wB*(1+half)
    grid = np.unique(np.concatenate([
        np.linspace(lo_w, hi_w, n),
        wB*(1 + np.concatenate([-np.logspace(-7, -1.2, 40), np.logspace(-7, -1.2, 40)])),
        [wB]]))
    P = np.array([prop(path, w, wr, mr) for w in grid])
    if not P.any(): return None                      # w_B sits inside one big stop band
    i0 = int(np.argmin(np.abs(grid - wB)))
    if not P[i0]:                                    # w_B not propagating: take the nearest run
        cand = np.where(P)[0]
        i0 = cand[int(np.argmin(np.abs(cand - i0)))]
    a = i0
    while a > 0 and P[a-1]: a -= 1
    b = i0
    while b < len(P)-1 and P[b+1]: b += 1
    if a == 0 or b == len(P)-1: return ('unbounded', grid[a], grid[b])
    def refine(x0, x1, want):
        while x1 - x0 > tol:
            m = 0.5*(x0+x1)
            if prop(path, m, wr, mr) == want: x0 = m
            else: x1 = m
        return 0.5*(x0+x1)
    lo = refine(grid[a], grid[a-1], True) if grid[a] > grid[a-1] else refine(grid[a-1], grid[a], False)
    lo = refine(grid[a-1], grid[a], False)
    hi = refine(grid[b], grid[b+1], True)
    return ('bounded', lo, hi)

if __name__ == "__main__":
    import time
    for kind, wB in (('membrane', np.pi), ('plate', np.pi**2)):
        pa = Path(kind)
        t0 = time.time(); prop(pa, wB*0.9, wB*0.8, 1.0)
        print(f"\n{'='*76}\n{kind.upper()}  w_B = {wB:.6f}   "
              f"[{(time.time()-t0)*1000:.1f} ms per propagation test]")
        for mr in (1.0, 4.0):
            print(f"   m_r = {mr}")
            for ratio in (0.90, 0.97, 1.00, 1.03, 1.10):
                wr = ratio*wB
                r = passband_around(pa, wB, wr, mr)
                bel, abv = edge_or_interior(pa, wB, wr, mr)
                cls = ("interior to a gap" if not any(bel+abv) else
                       "edge" if any(bel) != any(abv) else "interior to a passband")
                if r is None:
                    print(f"      w_r/w_B={ratio:.2f}  no propagating state near w_B   [{cls}]")
                elif r[0] == 'unbounded':
                    print(f"      w_r/w_B={ratio:.2f}  passband reaches the window edge "
                          f"({r[1]:.4f}-{r[2]:.4f})   [{cls}]")
                else:
                    print(f"      w_r/w_B={ratio:.2f}  passband {r[1]:9.5f} - {r[2]:9.5f}"
                          f"   width {r[2]-r[1]:8.5f} ({(r[2]-r[1])/wB*100:5.2f}% w_B)   [{cls}]")
