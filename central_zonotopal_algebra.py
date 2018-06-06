from abstract_zonotopal_algebra import AbstractZonotopalAlgebra

import poly_utils
from poly_free_module import PolynomialFreeModule

from sage.misc.cachefunc import cached_method

import itertools

# TODO Refactor functions to give dictionaries identifying polynomials with
#      corresponding matroid subsets, e.g. hyperplanes or bases, indep. sets


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
        basis = {}
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        X_cols = self.matrix().columns()
        for B in M.bases():
            ext_passive = M.passive_elements(B) - B
            elt = poly_utils.pure_tensor(P, X_cols, ext_passive)
            basis[B] = elt
        return basis

    @cached_method
    def D_space_basis(self):
        basis = {}
        P = self.polynomial_ring()
        M = self._ordered_matroid()
        ord_groundset = M.groundset_order()
        # handle ordering convention in OrderedMatroid class
        ord_groundset.reverse()

        def gs_key(x):
            return M.size() - M._gs_key(x)

        ext_ord = M.external_order(variant='convex geometry',
                                   representation='independent')
        X_cols = self.matrix().columns()

        # Cache dominant bases of each flat
        # empty flat
        dom_bases = {}
        dom_bases[frozenset([])] = frozenset([])
        # nonempty flats
        nonempty_flats = itertools.chain(
            *[M.flats(i) for i in range(1, M.rank() + 1)])
        for F in nonempty_flats:
            dom_bases[F] = M._dominant_basis(F)

        # For each independent set I, construct the D-space basis polynomial by
        # extending the basis polynomial associated with I - x where x is the
        # maximal element of I

        # base case: empty set
        basis[frozenset([])] = P.one()

        # recursively construct for additional elements in ord_groundset
        for x in ord_groundset:
            basis_update = {}
            for I0 in basis:
                # TODO give a more efficient enumeration of ind. set extensions
                I = I0 | frozenset([x])
                if not M.is_independent(I):
                    continue

                # identify starting D-space polynomial
                d0 = basis[I0]

                # identify flat for computation and J-generator for differentiation
                F0 = M.closure(I0)
                F = M.closure(I)
                # note for comparisons that the reverse order is used for notation
                cocirc = frozenset(
                    filter(lambda c: gs_key(c) <= gs_key(x), F - F0)
                )
                J_gen = poly_utils.pure_tensor(P, X_cols, cocirc)

                # compute orthogonal polynomial p_eta
                orthog_vec = self._hyperplane_normal(F0, F)
                p_eta = poly_utils.linear_form(P, orthog_vec)

                # extend d0 by power of orthogonal vector
                d = d0 * p_eta**(J_gen.degree() - 1)

                # compute derivative of d1 by J_gen
                d_deriv = poly_utils.poly_deriv(J_gen, d)

                # only project if this derivative is nonzero
                if d_deriv == P.zero():
                    basis_update[I] = d
                else:
                    # construct polynomial vector space for projection
                    dom_basis = dom_bases[F0]
                    poly_indices = ext_ord.closed_interval(I0, dom_basis)
                    poly_indices.remove(I0)
                    polys = [p_eta**(d.degree() - basis[J].degree()) * basis[J]
                             for J in poly_indices]
                    poly_derivs = [poly_utils.poly_deriv(J_gen, p)
                                   for p in polys]
                    P_mod = PolynomialFreeModule(P, basis=tuple(poly_derivs))

                    # decompose d derivative in this polynomial vector space
                    decomposition = P_mod(d_deriv).to_vector()
                    d_proj = d
                    for coeff, poly in zip(decomposition, polys):
                        d_proj -= coeff * poly
                    basis_update[I] = d_proj

            # update basis with new polynomials including x
            basis.update(basis_update)

        # normalize polynomials against P-space
        D_basis = {}
        P_basis = self.P_space_basis()
        for B in M.bases():
            coeff = poly_utils.diff_bilinear_form(P_basis[B], basis[B])
            D_basis[B] = basis[B] / coeff

        return D_basis
