from ordered_matroid import OrderedMatroid
from poly_free_module import PolynomialFreeModule

from sage.matroids.constructor import Matroid
from sage.modules.free_module import VectorSpace
from sage.rings.ideal import Ideal
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class AbstractZonotopalAlgebra:
    def __init__(self, X, varNames="x"):
        self._F = X.base_ring()
        self._M = Matroid(matrix=X)
        self._OM = OrderedMatroid(self._M, reverse=True)
        self._X = self._M.representation()
        self._V = VectorSpace(self._F, self._X.nrows())
        self._Pi = PolynomialRing(
            self._F, self._X.nrows(), names=varNames, order='deglex')

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

    def _ordered_matroid(self):
        return self._OM

    def _hyperplane_normal(self, H):
        X, V = self.matrix(), self._vector_space()
        H_subspace = V.subspace([X.column(h) for h in H])
        return H_subspace.complement().basis()[0]

    # TODO implement Lenz's generalized forward exchange cocircuits

    def I_ideal(self):
        return Ideal(self.I_ideal_gens())

    def J_ideal(self):
        return Ideal(self.J_ideal_gens())

    def I_ideal_gens(self):
        raise NotImplementedError

    def J_ideal_gens(self):
        raise NotImplementedError

    def P_space(self):
        return PolynomialFreeModule(self.polynomial_ring(),
                                    basis=tuple(self.P_space_basis()))

    def P_space_basis(self):
        raise NotImplementedError

    def D_space(self):
        raise PolynomialFreeModule(self.polynomial_ring(),
                                   basis=tuple(self.D_space_basis()))

    def D_space_basis(self):
        raise NotImplementedError
