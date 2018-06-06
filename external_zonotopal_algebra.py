from abstract_zonotopal_algebra import AbstractZonotopalAlgebra
from central_zonotopal_algebra import CentralZonotopalAlgebra

import poly_utils

from sage.matrix.constructor import Matrix
from sage.misc.cachefunc import cached_method


class ExternalZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self, X, varNames="x", externalBasisMatrix=None):
        if externalBasisMatrix is None:
            externalBasisMatrix = X.column_space().basis_matrix().transpose()
        else:
            if X.column_space() != externalBasisMatrix.column_space():
                raise ValueError(
                    "ExternalZonotopalAlgebra: externalBasisMatrix must have"
                    " the same column space as X")

        AbstractZonotopalAlgebra.__init__(self, X, varNames)
        self._ext_basis_matrix = externalBasisMatrix
        self._ext_block_matrix = Matrix.block([[X, self._ext_basis_matrix]])
        self._embedding_central_za = CentralZonotopalAlgebra(
            self._ext_block_matrix, varNames)

    def __repr__(self):
        return "External Zonotopal Algebra over " + str(self.base_field()) \
            + " with matrix\n" + str(self.matrix())

    def external_matrix(self):
        return self._ext_block_matrix

    def _external_matroid(self):
        return self._embedding_central_za._ordered_matroid()

    def _extending_basis_matrix(self):
        return self._ext_basis_matrix

    def _extending_basis_elements(self):
        r"""
        Return the extending basis elements as matroid groundset elements of
        the external matroid.
        """
        G = self._matroid().groundset()
        G_ext = self._external_matroid().groundset()
        return G_ext - G

    def _external_basis(self, I):
        r"""
        Return the external basis associated with a given independent set.
        """
        M_ext = self._external_matroid()
        rank = self.matrix().rank()
        extending_basis = sorted(self._extending_basis_elements())
        B_ext = tuple(I)
        # greedily add elements from extending_basis until full rank
        for b in extending_basis:
            if len(B_ext) == rank:
                break
            B_try = B_ext + (b,)
            if M_ext.is_independent(B_try):
                B_ext = B_try
        return frozenset(B_ext)

    def _external_bases(self):
        M = self._matroid()
        for I in M.independent_sets():
            yield self._external_basis(I)

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

    @cached_method
    def J_ideal_gens(self):
        P = self.polynomial_ring()
        M = self._external_matroid()
        X_cols = self.external_matrix().columns()
        # compute internal generalized cocircuits
        L = M.external_order(variant="antimatroid")
        non_external = set(L) - set(self._external_bases())
        min_non_external = L.subposet(non_external).minimal_elements()
        gen_cocircuits = [M.passive_elements(I) - I for I in min_non_external]
        # initialize generators
        gens = []
        for c in gen_cocircuits:
            gen = poly_utils.pure_tensor(P, X_cols, c)
            gens.append(gen)
        return gens

    @cached_method
    def P_space_basis(self):
        basis = {}
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        X_cols = self.matrix().columns()
        for I in M.independent_sets():
            ext_passive = M.passive_elements(I) - I
            elt = poly_utils.pure_tensor(P, X_cols, ext_passive)
            basis[I] = elt
        return basis

    @cached_method
    def D_space_basis(self):
        M = self._matroid()
        central_basis = self._embedding_central_za.D_space_basis()
        basis = {I: central_basis[self._external_basis(I)]
                 for I in M.independent_sets()}
        return basis
