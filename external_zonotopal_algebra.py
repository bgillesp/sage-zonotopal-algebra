from abstract_zonotopal_algebra import AbstractZonotopalAlgebra
from poly_utils import PolyUtils
from sage.misc.cachefunc import cached_method
from sage.matroids.constructor import Matroid
# from sage.modules.free_module import VectorSpace
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.misc.misc_c import prod
# from sage.functions.other import factorial
from sage.matrix.constructor import Matrix
from matroid_utils import cocircuits

class ExternalZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self,X,varNames="x",externalBasisMatrix=None):
        if externalBasisMatrix == None:
            externalBasisMatrix = Matrix.identity(X.nrows())
        AbstractZonotopalAlgebra.__init__(self, X, varNames)
        self._ext_basis_matrix = externalBasisMatrix
        self._ext_block_matrix = Matrix.block([[X, self._ext_basis_matrix]])
        self._ext_matroid = Matroid(self._ext_block_matrix)

    def __repr__(self):
        return "External Zonotopal Algebra over " + str(self.base_field()) \
            + " with matrix\n" + str(self.matrix())

    # TODO is it better to keep track of base matroid, or matroid with extra basis?

    def external_matrix(self):
        return self._ext_block_matrix

    def _external_basis(self):
        return self._ext_basis_matrix

    def _external_matroid(self):
        return self._ext_matroid

    def _external_bases(self):
        M0 = self._external_matroid()
        rank = self.matrix().nrows()
        extendors = list(M0.groundset())[-rank:]
        for x in self._matroid().independent_sets():
            x = list(x)
            for y in extendors:
                if len(x) == rank:
                    break
                if M0.is_independent(x + [y]):
                    x = x + [y]
            yield x

    @cached_method
    def I_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        n = self._matroid().size()
        for normal,hyperplane in zip( self._facet_hyperplane_normals(), self._facet_hyperplanes() ):
            gen = PolyUtils.linear_form(P, normal)**(n - len(hyperplane) + 1)
            gens.append(gen)
        return gens

    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X_cols = self.external_matrix().columns()
        for c in cocircuits(self._external_matroid(), self._external_bases()):
            gen = PolyUtils.pure_tensor(P, X_cols, c)
            gens.append(gen)
        return gens

    @cached_method
    def P_basis(self):
        basis = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        for indep in self._matroid().independent_sets():
            x_indep = self._big_ex(indep)
            elt = PolyUtils.pure_tensor(P, X_cols, x_indep)
            basis.append(elt)
        return basis

    # TODO implement D_basis
    # @cached_method
    # def D_basis(self):
    #     basis = []
    #     P = self.polynomial_ring()
    #     X_cols = self.matrix().columns()
    #     # compute dual of P(X)
    #     P_dual_basis = PolyUtils.poly_dual_basis(P, self.P_basis())
    #     # project P(X) basis (acting as a P(X)* basis) into the kernel of J(X)
    #     # (hyperplane normals and cocircuits are ordered to allow zipping)
    #     polys = []
    #     for normal, cocirc in zip(self._facet_hyperplane_normals(), self._cocircuits()):
    #         i = PolyUtils.linear_form(P, normal)**len(cocirc)
    #         j = PolyUtils.pure_tensor(P, X_cols, cocirc)
    #         polys.append( (i,j) )
    #     for p in P_dual_basis:
    #         q = p
    #         for (i,j) in polys:
    #             coeff1 = PolyUtils.diff_bilinear_form(j,p)
    #             coeff2 = PolyUtils.diff_bilinear_form(j,i)
    #             q -= coeff1/coeff2 * i
    #         basis.append(q)
    #     return basis
