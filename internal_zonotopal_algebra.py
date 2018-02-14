# from polynomial_vector_space import PolynomialModule
# from polynomial_vector_space import Monomials
# from sage.misc.cachefunc import cached_method
# from sage.matroids.constructor import Matroid
# from sage.modules.free_module import VectorSpace
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.misc.misc_c import prod
# from sage.functions.other import factorial
# from sage.matrix.constructor import Matrix
import poly_utils
# from set_utils import SetUtils
import matroid_utils
from sage.misc.cachefunc import cached_method
from abstract_zonotopal_algebra import AbstractZonotopalAlgebra


class InternalZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self, X, varNames="x"):
        AbstractZonotopalAlgebra.__init__(self, X, varNames)

    def __repr__(self):
        return ("Internal Zonotopal Algebra over "
                + str(self.base_field())
                + " with matrix\n"
                + str(self.matrix()))

    def _internal_bases(self):
        M = self._ordered_matroid()
        G = M.groundset()
        for B in M.bases():
            int_active = M.coactive_elements(G - B) & B
            # no internally active elements
            if len(int_active) == 0:
                yield B

    @cached_method
    def I_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        n = self._matroid().size()
        for H in self._matroid().hyperplanes():
            normal = self._hyperplane_normal(H)
            gen = poly_utils.linear_form(P, normal)**(n - len(H) - 1)
            gens.append(gen)
        return gens

    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        for c in matroid_utils.cocircuits(self._matroid(),
                                          self._internal_bases()):
            gen = poly_utils.pure_tensor(P, X_cols, c)
            gens.append(gen)
        return gens

    @cached_method
    def P_basis(self):
        basis = []
        # for each element of internally passive bases, check if ext active set
        # in cocircuit is empty
        # if so, zero out b-component of largest elt in ext passive set
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        X = self.matrix()
        X_cols = self.matrix().columns()
        for B in self._internal_bases():
            ext_passive = M.passive_elements(B) - B
            projections = []
            for b in B:
                fund_cocirc = M.fundamental_cocircuit(B, b)
                if len(fund_cocirc.difference(ext_passive)) == 1:
                    # only b is there, no externally active elts
                    # project maximal ext passiv elt in X_b away from b
                    projected = max(fund_cocirc.intersection(ext_passive))
                    projections.append((b, projected))

            if len(projections) > 0:
                # 1. represent vectors in ext_passive in terms of basis B
                basis_matrix = X[:, list(B)]
                basis_indices = {b: i for i, b in enumerate(B)}
                X_basis_change = basis_matrix.solve_right(X)
                # 2. project passive vectors away from corresponding basis elts
                for b, proj in projections:
                    X_basis_change[basis_indices[b], proj] = 0
                # 3. represent projected matrix in terms of ambient basis again
                X_projected = basis_matrix * X_basis_change

                elt = poly_utils.pure_tensor(
                    P, X_projected.columns(), ext_passive)
            else:
                elt = poly_utils.pure_tensor(P, X_cols, ext_passive)
            basis.append(elt)
        return basis
