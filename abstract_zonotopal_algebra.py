from poly_utils import PolyUtils
# from sage.misc.cachefunc import cached_method
from sage.matroids.constructor import Matroid
from sage.modules.free_module import VectorSpace
from sage.rings.ideal import Ideal
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod

class AbstractZonotopalAlgebra:
    def __init__(self, X, varNames="x"):
        self._F = X.base_ring()
        self._M = Matroid(matrix=X)
        self._X = self._M.representation()
        self._V = VectorSpace(self._F,self._X.nrows())
        self._Pi = PolynomialRing(self._F, self._X.nrows(), names=varNames, order='deglex')

    def base_field(self):
        return self._F

    def polynomial_ring(self):
        return self._Pi

    def matrix(self):
        return self._X

    def _vector_space(self):
        return self._V

    def _matroid(self):
        return self._M

    def _facet_hyperplanes(self):
        # Matroid.hyperplanes() should be a (lazy) iterator for efficiency
        for h in self._matroid().hyperplanes():
            yield h

    def _facet_hyperplane_normals(self):
        X = self.matrix()
        M = self._matroid()
        V = self._vector_space()
        for h in self._facet_hyperplanes():
            h_basis = [X.column(v) for v in M.max_independent(h)]
            yield V.subspace(h_basis).complement().basis()[0]

    def _cocircuits(self):
        g = self._matroid().groundset()
        for h in self._facet_hyperplanes():
            yield g.difference(h)

    # TODO implement generalized cocircuits for forward exchange matroids, from Lenz

    # TODO is this abstract, or specific to central case?
    def _big_ex(self,S):
        M = self._matroid()
        ex_poss = M.groundset().difference(S)
        def f(e):
            flat = M.closure(filter(lambda i: i <= e, S))
            return not e in flat
        ex = filter(f, ex_poss)
        return ex

    # TODO what is this?
    def _big_y(self,s):
        M = self._matroid()
        y_poss = M.groundset().difference(s)
        def f(e):
            flat = M.closure(filter(lambda i: i >= e, s))
            return not e in flat
        y = filter(f, y_poss)
        return y

    def I_ideal(self):
        return Ideal(self.I_ideal_gens())

    def J_ideal(self):
        return Ideal(self.J_ideal_gens())

    def I_ideal_gens(self):
        raise NotImplementedError

    def J_ideal_gens(self):
        raise NotImplementedError

    def P_basis(self):
        raise NotImplementedError

    def D_basis(self):
        raise NotImplementedError
