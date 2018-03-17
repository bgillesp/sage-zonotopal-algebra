from abstract_zonotopal_algebra import AbstractZonotopalAlgebra

import poly_utils

from sage.misc.cachefunc import cached_method


class CentralZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self, X, varNames="x"):
        AbstractZonotopalAlgebra.__init__(self, X, varNames)

    def __repr__(self):
        return "Central Zonotopal Algebra over " + str(self.base_field()) \
            + " with matrix\n" + str(self.matrix())

    @cached_method
    def I_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        n = self._matroid().size()
        for H in self._matroid().hyperplanes():
            normal = self._hyperplane_normal(H)
            gen = poly_utils.linear_form(P, normal)**(n - len(H))
            gens.append(gen)
        return gens

    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        for cocirc in self._matroid().cocircuits():
            gen = poly_utils.pure_tensor(P, X_cols, cocirc)
            gens.append(gen)
        return gens

    @cached_method
    def P_space_basis(self):
        basis = []
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        X_cols = self.matrix().columns()
        for b in M.bases():
            ext_passive = M.passive_elements(b) - b
            elt = poly_utils.pure_tensor(P, X_cols, ext_passive)
            basis.append(elt)
        return basis

    @cached_method
    def D_space_basis(self):
        basis = []
        P = self.polynomial_ring()
        M = self._matroid()
        X_cols = self.matrix().columns()
        # compute dual of P(X) basis inside of P(X)
        P_dual_basis = poly_utils.poly_dual_basis(P, self.P_basis())
        # project dual basis into the kernel of J(X)
        polys = []
        for H in M.hyperplanes():
            normal = self._hyperplane_normal(H)
            cocirc = M.groundset() - H
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
