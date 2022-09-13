from typing import Optional, List
from sage.all import *
import math
import numpy as np

class SELCSolver:
    def __init__(self, maxn: int, deg: Optional[int] = None, iter: Optional[int] = None) -> None:
        
        self.maxn: int = maxn

        self.deg: int = 0
        if deg is None:
            self.deg = math.ceil(5 * math.log2(maxn))
        else:
            self.deg: int = deg

        self.iter: int = 0
        if iter is None:
            self.iter = self.maxn + 1
        else:
            self.iter = iter

        # finite field
        self.F = GF(Integer(2 ** self.deg), names=('x',))
        # print(f"Field degree: {self.F.degree()}")
        (self.x,) = self.F._first_ngens(1)

        self.Fx = self.F['t']
        (self.t,) = self.Fx._first_ngens(1)

        # extension ring
        self.E = PolynomialRing(Integers(Integer(4)), 'x').quotient(self.F.modulus(), names=('y',))
        (self.y,) = self.E._first_ngens(1)

        self.perm_call_count = 0

    def source_distinct_elements_from_field(self, n: int) -> List:
        S = list(set([self.F.random_element() for _ in range(2 * n)]))
        assert len(S) >= n, f"failed to source distinct {n} elements from GF(2^{self.deg})"
        return S[:n]

    '''
    returns an instance of the shortest even cycle, or [] if it does not exist
    '''
    def get_selc(self, a_: np.ndarray) -> List[int]:
        self.perm_call_count = 0
        a = a_.copy().astype(np.int32)
        n = a.shape[0]
        selc = self.get_selc_length(a)
        # print(f"selc length: {selc}")

        if selc == -1:
            return []
        
        ver = list(range(n))

        for i in range(n-1, -1, -1):
            if len(ver) == selc:
                break

            b = np.delete(np.delete(a, i, axis=1), i, axis=0)
            if self.get_selc_length(b) == selc:
                ver.remove(i)
                a = b

        cycle = []
        now = 0
        while len(cycle) < len(ver) - 1:
            (out,) = np.nonzero(a[now])
            s, e = 0, len(out) - 1
            
            while s < e:
                m = (s + e) >> 1
                filt = np.zeros(len(ver)).astype(np.int32)
                filt[out[s] : out[m] + 1] = 1
                a[now] *= filt

                if self.get_selc_length(a) == selc:
                    e = m
                else:
                    s = m + 1
            
            cycle.append(ver[now])
            now = out[s]

        cycle.append(ver[now])
        return cycle
            

    '''
    returns length of the shortest even cycle if exists, otherwise -1.
    '''
    def get_selc_length(self, a: np.ndarray) -> int:
        assert a.ndim == 2, "only 2-dimensional matrices are allowed"
        assert a.shape[0] == a.shape[1], "square adjacency matrix is required"
        assert a.dtype in [np.int32, np.int64], "integral adjacency matrix is required"
        n = a.shape[0]

        ans = n + 1
        gamma = self.source_distinct_elements_from_field(n + 1)

        for _ in range(self.iter):
            ans = min(ans, self.oracle(a, gamma))
        
        return -1 if ans > n else ans

    # def get_odd_elim_multiplier(self, src, dst):
    #     print(src, dst, src.list())
    #     fs = self.F(src.list())
    #     assert fs != 0, "src must be odd"
    #     return self.E(1 / fs) * dst

    def permanent_with_similar_rows_lagrange(self, mat, i1, i2):
        self.perm_call_count += 1
        def projection(e):
            return self.F(list(map( lambda x: x % 2, e.list() )))

        n = mat.nrows()
        M = Matrix(self.F,
        [[ projection(mat[i, j]) for j in range(n) ] for i in range(n)]
        )

        points = self.source_distinct_elements_from_field(2*n - 1)
        
        itp_witness = []
        for x in points:
            N = copy(M)
            for j in range(n):
                N[i1, j] *= x ** Integer(j)
                N[i2, j] *= x ** Integer(n-1-j)
            
            itp_witness.append((x, N.det()))
        
        pol = self.Fx.lagrange_polynomial(itp_witness)
        return self.E(sum(pol.list()[:n-1])) * 2

    def permanent_with_similar_rows(self, mat, i1, i2):

        def projection(e):
            return self.F(list(map( lambda x: x % 2, e.list() )))

        n = mat.nrows()
        M = Matrix(self.Fx,
        [[ projection(mat[i, j]) for j in range(n) ] for i in range(n)]
        )
        for j in range(n):
            M[i1, j] *= self.t ** Integer(j)
            M[i2, j] *= self.t ** Integer(n-1-j)
        
        return self.E(sum(M.det().list()[:n-1])) * 2


    def perm_bjorklund(self, mat):
        M = copy(mat)

        n = mat.nrows()
        max_marked_row = -1
        max_marked_col = -1
        def get_col_with_odd_unmarked_row():
            for j in range(max_marked_col + 1, n):
                if (2 * M[max_marked_row + 1 :, j]).is_zero():
                    continue
                for i in range(max_marked_row + 1, n):
                    if 2 * M[i, j] != 0:
                        return i, j
            return None, None

        tot_perm = self.E(0)
        while True:
            i1, j = get_col_with_odd_unmarked_row()
            if i1 is None:
                break

            inv_of_odd = self.E(self.F.one() / self.F(list(map(lambda x : x % 2, M[i1, j].list()))))
            # print(inv_of_odd)
            for i2 in range(n):
                if i1 == i2:
                    continue

                if 2 * M[i2, j] != 0:
                    rho = inv_of_odd * M[i2, j]
                    N = copy(M)
                    N.rescale_row(i2, 0)
                    N.add_multiple_of_row(i2, i1, rho)

                    # from time import perf_counter as pf
                    # t = pf()
                    # perm1 = self.permanent_with_similar_rows(N, i1, i2)
                    # t2= pf()
                    perm2 = self.permanent_with_similar_rows_lagrange(N, i1, i2)
                    # t3 = pf()
                    # assert perm1 == perm2, f"perm1 = {perm1}, perm2 = {perm2}"
                    # print(f"naive: {t2 - t:.3f}, ho: {t3 - t2:.3f}")
                    tot_perm += perm2
                    M.add_multiple_of_row(i2, i1, rho)
            
            M.swap_rows(i1, max_marked_row + 1)
            M.swap_columns(j, max_marked_col + 1)
            max_marked_col += 1
            max_marked_row += 1
        
        if max_marked_row < n-2:
            return tot_perm
        
        return tot_perm + product(M.diagonal())


    def perm_sage(self, mat):
        return mat.permanent()

    def det_sage(self, mat):
        return mat.det()
    
    def oracle(self, a: np.ndarray, gamma: List) -> int:
        n = a.shape[0]

        interpolation_points = []        
        beta = random_matrix(self.F, n, n)

        for gam in gamma:
            a_in_f = Matrix(self.F, a.tolist())
            for i in range(n):
                a_in_f[i, i] = gam

            mat = Matrix(self.E, beta.elementwise_product(a_in_f))

            perm = self.perm_bjorklund(mat)
            pcc_double = perm - self.det_sage(mat)

            # psage = self.perm_sage(mat)
            # assert perm == psage, f"Poly: {perm}, Sage: {psage}"
            
            pcc_coeffs = pcc_double.list()
            for coeff in pcc_coeffs:
                if coeff % 2:
                    assert False, f"{pcc_coeffs} has an odd coefficient"
            
            pcc_coeffs = list(map(lambda x : int(x) // 2, pcc_coeffs))
            pcc = self.F(pcc_coeffs)

            interpolation_points.append((gam, pcc))
        
        Q = self.Fx.lagrange_polynomial(interpolation_points)
        # print(Q)

        ql = Q.list()
        for k in range(2, n+1, 2):
            if n - k < len(ql) and ql[n-k] != 0:
                # print(f"SELC is {k}")
                return k
        
        return n + 1