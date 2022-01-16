from sage.matroids.constructor import Matroid
from sage.modules.free_module import VectorSpace
from sage.rings.ideal import Ideal
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from .ordered_matroid import OrderedMatroid
from .poly_free_module import PolynomialFreeModule


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

    def _hyperplane_normal(self, F, E=None):
        r"""
        Return a hyperplane normal of a flat F in another flat E of rank one
        higher.

        INPUT:

        - ``F`` -- a flat of the underlying matroid.

        - ``E`` -- (default: ``None``) a flat of the underlying matroid of rank
          one higher than ``F``.  If ``None``, then ``E`` is assumed to be the
          complete ground set.

        OUTPUT:

        A nonzero vector in the span of ``E`` which is orthogonal to the span
        of ``F``.
        """
        X, V = self.matrix(), self._vector_space()

        F_subspace = V.subspace([X.column(x) for x in F])

        if E is None:
            E = self._matroid().groundset()
            # return F_subspace.complement().basis()[0]
        E_subspace = V.subspace([X.column(x) for x in E])

        # Workaround needed for bug similar to the one reported at:
        # https://trac.sagemath.org/ticket/32447
        # Supposedly should be fixed in Sage 9.5
        #
        # Original code:
        # return F_subspace.complement().intersection(E_subspace).basis()[0]

        external_vec = None
        for x in E:
            if X.column(x) not in F_subspace:
                external_vec = X.column(x)
                break
        basis = F_subspace.basis()
        for b in basis:
            external_vec -= \
                (b.inner_product(external_vec) / b.inner_product(b)) * b

        return external_vec



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
