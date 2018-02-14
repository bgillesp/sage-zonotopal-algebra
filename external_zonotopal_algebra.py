from abstract_zonotopal_algebra import AbstractZonotopalAlgebra
import poly_utils
from sage.misc.cachefunc import cached_method
from sage.matroids.constructor import Matroid
# from sage.modules.free_module import VectorSpace
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.misc.misc_c import prod
# from sage.functions.other import factorial
from sage.matrix.constructor import Matrix
from matroid_utils import cocircuits


class ExternalZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self, X, varNames="x", externalBasisMatrix=None):
        if externalBasisMatrix is None:
            externalBasisMatrix = Matrix.identity(X.nrows())
        AbstractZonotopalAlgebra.__init__(self, X, varNames)
        self._ext_basis_matrix = externalBasisMatrix
        self._ext_block_matrix = Matrix.block([[X, self._ext_basis_matrix]])
        self._ext_matroid = Matroid(self._ext_block_matrix)

    def __repr__(self):
        return "External Zonotopal Algebra over " + str(self.base_field()) \
            + " with matrix\n" + str(self.matrix())

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
        for H in self._matroid().hyperplanes():
            normal = self._hyperplane_normal(H)
            gen = poly_utils.linear_form(P, normal)**(n - len(H) + 1)
            gens.append(gen)
        return gens

    # def alt_I_ideal_gens(self):
    #     # IDEA: generate I ideal using some appropriate generators coming from
    #     # lower-dimensional subspaces, not just the hyperplanes
    #     # Take k normals, forming a flag in the orthogonal space
    #     # Take a product of powers of the normal linear forms corresponding to
    #     # vectors outside of the subspaces spanned by successively more of the normals
    #     # I think the theory works out for this, since you'll get a nonzero bilinear form for the
    #     # corresponding cocircuit generator of the external J-ideal.
    #     # - Are these elements also contained in the standard external I-ideal?
    #     # - Where does the +1 come in in the exponent?
    #     gens = []
    #     # normals for arbitrary subspaces
    #     X = self.matrix()
    #     M = self._matroid()
    #     G = M.groundset()
    #     V = self._vector_space()
    #     flat_bases = []
    #     for c in cocircuits(self._external_matroid(), self._external_bases()):
    #         flat = G.difference(c)
    #         f_basis = [X.column(v) for v in M.max_independent(flat)]
    #         f_orth_basis = V.subspace(f_basis).complement().basis()
    #         flat_bases.append( (flat, f_orth_basis) )
    #     return flat_bases

    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X_cols = self.external_matrix().columns()
        for c in cocircuits(self._external_matroid(), self._external_bases()):
            gen = poly_utils.pure_tensor(P, X_cols, c)
            gens.append(gen)
        return gens

    @cached_method
    def P_basis(self):
        basis = []
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        X_cols = self.matrix().columns()
        for indep in M.independent_sets():
            ext_passive = M.passive_elements(indep) - indep
            elt = poly_utils.pure_tensor(P, X_cols, ext_passive)
            basis.append(elt)
        return basis

    # TODO implement D_basis
    # @cached_method
    # def D_basis(self):
    #     basis = []
    #     P = self.polynomial_ring()
    #     X_cols = self.matrix().columns()
    #     # compute dual of P(X)
    #     P_dual_basis = poly_utils.poly_dual_basis(P, self.P_basis())
    #     # project P(X) basis (acting as a P(X)* basis) into the kernel of J(X)
    #     # (hyperplane normals and cocircuits are ordered to allow zipping)
    #     polys = []
    #     for normal, cocirc in zip(self._facet_hyperplane_normals(), self._cocircuits()):
    #         i = poly_utils.linear_form(P, normal)**len(cocirc)
    #         j = poly_utils.pure_tensor(P, X_cols, cocirc)
    #         polys.append( (i,j) )
    #     for p in P_dual_basis:
    #         q = p
    #         for (i,j) in polys:
    #             coeff1 = poly_utils.diff_bilinear_form(j,p)
    #             coeff2 = poly_utils.diff_bilinear_form(j,i)
    #             q -= coeff1/coeff2 * i
    #         basis.append(q)
    #     return basis
