from abstract_zonotopal_algebra import AbstractZonotopalAlgebra
from poly_utils import PolyUtils
from sage.misc.cachefunc import cached_method
# from sage.matroids.constructor import Matroid
# from sage.modules.free_module import VectorSpace
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.misc.misc_c import prod
# from sage.functions.other import factorial
# from sage.matrix.constructor import Matrix

class CentralZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self,X,varNames="x"):
        AbstractZonotopalAlgebra.__init__(self, X, varNames)

    def __repr__(self):
        return "Central Zonotopal Algebra over " + str(self.base_field()) \
            + " with matrix\n" + str(self.matrix())

    @cached_method
    def I_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        n = self._matroid().size()
        for normal,hyperplane in zip( self._facet_hyperplane_normals(), self._facet_hyperplanes() ):
            gen = PolyUtils.linear_form(P, normal)**(n - len(hyperplane))
            gens.append(gen)
        return gens

    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        for cocirc in self._cocircuits():
            gen = PolyUtils.pure_tensor(P, X_cols, cocirc)
            gens.append(gen)
        return gens

    @cached_method
    def P_basis(self):
        basis = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        for b in self._matroid().bases():
            xb = self._big_ex(b)
            elt = PolyUtils.pure_tensor(P, X_cols, xb)
            basis.append(elt)
        return basis

    @cached_method
    def D_basis(self):
        basis = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        # compute dual of P(X)
        P_dual_basis = PolyUtils.poly_dual_basis(P, self.P_basis())
        # project P(X) basis (acting as a P(X)* basis) into the kernel of J(X)
        # (hyperplane normals and cocircuits are ordered to allow zipping)
        polys = []
        for normal, cocirc in zip(self._facet_hyperplane_normals(), self._cocircuits()):
            i = PolyUtils.linear_form(P, normal)**len(cocirc)
            j = PolyUtils.pure_tensor(P, X_cols, cocirc)
            polys.append( (i,j) )
        for p in P_dual_basis:
            q = p
            for (i,j) in polys:
                coeff1 = PolyUtils.diff_bilinear_form(j,p)
                coeff2 = PolyUtils.diff_bilinear_form(j,i)
                q -= coeff1/coeff2 * i
            basis.append(q)
        return basis
