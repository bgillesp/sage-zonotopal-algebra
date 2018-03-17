from abstract_zonotopal_algebra import AbstractZonotopalAlgebra

import poly_utils

from sage.matroids.constructor import Matroid
from sage.matrix.constructor import Matrix
from sage.misc.cachefunc import cached_method


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

    def _external_matroid(self):
        return self._ext_matroid

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
        rank = self.matrix().nrows()
        extending_basis = tuple(self._extending_basis_elements())
        B_ext = tuple(I)
        # greedily add elements from extending_basis until full rank
        for b in extending_basis:
            if len(B_ext) == rank:
                break
            B_try = B_ext + (b,)
            if M_ext.is_independent(B_try):
                B_ext = B_try
        return frozenset(B_ext)

    def _external_cocircuit(self, H):
        r"""
        Return the external cocircuit associated with a given hyperplane.

        The external cocircuit associated with a given hyperplane is the set
        obtained by greedily adding a single element of the extending basis
        not spanned by the hyperplane to the original cocircuit.
        """
        M_ext = self._external_matroid()
        if M_ext.rank(H) != M_ext.rank() - 1:
            raise ValueError("The given set H is not a hyperplane")
        H_ext = M_ext.closure(H)
        extending_basis = self._extending_basis_elements()
        G = self._matroid().groundset()
        extendor = frozenset(sorted(extending_basis - H_ext)[:1])
        H = frozenset(H)
        return (G - H_ext) | extendor

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
        gens = []
        P = self.polynomial_ring()
        X_cols = self.external_matrix().columns()
        for H in self._matroid().hyperplanes():
            cocirc = self._external_cocircuit(H)
            gen = poly_utils.pure_tensor(P, X_cols, cocirc)
            gens.append(gen)
        return gens

    @cached_method
    def P_space_basis(self):
        basis = []
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        X_cols = self.matrix().columns()
        for indep in M.independent_sets():
            ext_passive = M.passive_elements(indep) - indep
            elt = poly_utils.pure_tensor(P, X_cols, ext_passive)
            basis.append(elt)
        return basis

    @cached_method
    def D_space_basis(self):
        basis = []
        P = self.polynomial_ring()
        M = self._matroid()
        X_cols = self.external_matrix().columns()
        # compute dual of P(X) basis inside of P(X)
        P_dual_basis = poly_utils.poly_dual_basis(P, self.P_basis())
        # project dual basis into the kernel of J(X)
        polys = []
        for H in M.hyperplanes():
            normal = self._hyperplane_normal(H)
            cocirc = self._external_cocircuit(H)
            i = poly_utils.linear_form(P, normal)**len(cocirc)  # I-ideal gen
            j = poly_utils.pure_tensor(P, X_cols, cocirc)       # J-ideal gen
            polys.append((i, j))
        for p in P_dual_basis:
            q = p
            for (i, j) in polys:
                coeff1 = poly_utils.diff_bilinear_form(j, p)
                coeff2 = poly_utils.diff_bilinear_form(j, i)
                q -= (coeff1 / coeff2) * i
            basis.append(q)
        return basis
