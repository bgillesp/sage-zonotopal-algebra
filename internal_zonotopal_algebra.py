# from polynomial_vector_space import PolynomialModule
# from polynomial_vector_space import Monomials
# from poly_utils import PolyUtils
# from sage.misc.cachefunc import cached_method
# from sage.matroids.constructor import Matroid
# from sage.modules.free_module import VectorSpace
# from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
# from sage.misc.misc_c import prod
# from sage.functions.other import factorial
# from sage.matrix.constructor import Matrix
from poly_utils import PolyUtils
from sage.misc.cachefunc import cached_method
from abstract_zonotopal_algebra import AbstractZonotopalAlgebra

class InternalZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self,X,varNames="x"):
        AbstractZonotopalAlgebra.__init__(self, X, varNames)

    def __repr__(self):
        return "Internal Zonotopal Algebra over the field " + str(self.F) + " with matrix\n" + str(self.X)

    def _internal_bases(self):
        M = self._matroid()
        order = range(M.size())
        order.reverse() # formulate activity parameter in terms of max, as in (Holtz, Ron)
        for b in M.bases():
            if len(matroid_internal(M, b, order)) == 0: # no internally active elements
                yield b

    @cached_method
    def _min_long_subsets(self):
        return ZonUtils.minimal_sets(self.long_subsets())

    @cached_method
    def I_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        n = self._matroid().size()
        for normal,hyperplane in zip( self._facet_hyperplane_normals(), self._facet_hyperplanes() ):
            gen = PolyUtils.linear_form(P, normal)**(n - len(hyperplane) - 1)
            gens.append(gen)
        return gens

    # TODO check correctness of this
    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X = self.matrix()
        for c in self._min_long_subsets():
            gens.append(PolyUtils.pure_tensor(P, X.columns(), c))
        return gens

    @cached_method
    def P_basis(self):
        basis = []
        # for each element of internally passive bases, check if ext active set in cocircuit is empty
        # if so, zero out b-component of largest elt in ext passive set
        P = self.polynomial_ring()
        M = self._matroid()
        G = self._matroid().groundset()
        X = self.matrix()
        X_cols = self.matrix().columns()
        V = self._vector_space()
        F = self.base_field()
        for B in self._internal_bases():
            xb = self._big_ex(B) # externally passive elts
            projections = []
            for b in B:
                fund_cocirc = M.fundamental_cocircuit(B, b)
                if len(fund_cocirc.difference(xb)) == 1:
                    # only b is there, no externally active elts
                    # project maximal ext passiv elt in X_b away from b
                    projected = max( fund_cocirc.intersection(xb) )
                    projections.append( (b, projected) )

            if len(projections) > 0:

                ########### VERBOSE VERSION for debugging
                # # apply projections to corresponding vectors in copy of X
                # print "X:\n%s" % X
                # print "Basis: %s" % B
                # print "xb: %s" % xb
                # print "Projections: %s" % projections
                # X_copy = copy(X)
                # # 1. represent vectors in xb in terms of basis B
                # basis_matrix = X_copy[:,list(B)]
                # print "Basis matrix:\n%s" % basis_matrix
                # basis_indices = {b: i for i,b in enumerate(B)}
                # # could truncate X to only necessary columns, but simplifies
                # # code to leave it as is; modify this if efficiency is an issue
                # X_basis_change = basis_matrix.solve_right(X_copy)
                # print "X after basis change:\n%s" % X_basis_change
                # # 2. project xb vectors away from corresponding basis elts b
                # for b, proj in projections:
                #     X_basis_change[basis_indices[b], proj] = 0
                # print "Projected matrix:\n%s" % X_basis_change
                # # 3. represent again in terms of ambient basis
                # X_copy = basis_matrix * X_basis_change
                # print "X in ambient basis after projection:\n%s" % X_copy
                # # 4. translate to polynomial
                # elt = PolyUtils.pure_tensor(P, X_copy.columns(), xb)

                # apply projections to corresponding vectors in copy of X
                # 1. represent vectors in xb in terms of basis B
                basis_matrix = X[:,list(B)]
                basis_indices = {b: i for i,b in enumerate(B)}
                X_basis_change = basis_matrix.solve_right(X)
                # 2. project xb vectors away from corresponding basis elts b
                for b, proj in projections:
                    X_basis_change[basis_indices[b], proj] = 0
                # 3. represent projected matrix in terms of ambient basis again
                X_projected = basis_matrix * X_basis_change

                elt = PolyUtils.pure_tensor(P, X_projected.columns(), xb)
            else:
                elt = PolyUtils.pure_tensor(P, X_cols, xb)
            basis.append(elt)
        return basis

    # @cached_method
    # def D_basis(self):
    #     basis = []
    #     P = self.polynomial_ring()
    #     X = self.matrix()
    #     # compute dual of P(X)
    #     P_dual_basis = PolyUtils.poly_dual_basis(P, self.P_basis())
    #     # project P(X) basis (acting as a P(X)* basis) into the kernel of J(X)
    #     # (hyperplane normals and cocircuits are ordered to allow zipping)
    #     polys = []
    #     for normal, cocirc in zip(self._facet_hyperplane_normals(), self._cocircuits()):
    #         i = PolyUtils.linear_form(P, normal)**len(cocirc)
    #         j = PolyUtils.pure_tensor(P, X.columns, cocirc)
    #         polys.append( (i,j) )
    #     for p in P_dual_basis:
    #         q = p
    #         for (i,j) in polys:
    #             coeff1 = PolyUtils.diff_bilinear_form(j,p)
    #             coeff2 = PolyUtils.diff_bilinear_form(j,i)
    #             q -= coeff1/coeff2 * i
    #         basis.append(q)
    #     return basis


# TODO check these for correctness
def matroid_external(M, B, order=None):
    """
    Return the set of externally active elements of a basis `B`.

    An element `e` is *externally active* if it is the smallest element in
    the `B`-fundamental circuit using `e`. Smallest is interpreted as
    the output of the built-in ``min`` function on the subset.

    The `B`-fundamental circuit using `e` is the unique circuit contained
    in `B + e`.

    INPUT:

    - ``M`` -- an ambient matroid in which to compute the internally active
      elements
    - ``B`` -- a basis of the matroid, assumed to have Python's
      ``frozenset`` interface.
    - ``order`` -- an optional ordering of the base set of M to use
      for the computation.

    Function based on implementation in Matroid._internal.

    TODO:

    Might it be a good idea to make internal and external activity functions
    part of the public interface to the Matroid class?
    """
    N = M.groundset() - B
    A = set()
    order_key = lambda x: x
    if order != None:
        order_key = lambda x: order.index(x) # order should contain all indices
    for e in N:
        if min(M.circuit(B | set([e])), key=order_key) == e:
            A.add(e)
    return A

def matroid_internal(M, B, order=None):
    """
    Return the set of internally active elements of a basis `B`.

    An element `e` is *internally active* if it is the smallest element in
    the `B`-fundamental cocircuit using `e`. Smallest is interpreted as
    the output of the built-in ``min`` function on the subset.

    The `B`-fundamental cocircuit using `e` is the unique cocircuit
    intersecting basis `B` in exactly element `e`.

    INPUT:

    - ``M`` -- an ambient matroid in which to compute the internally active
      elements
    - ``B`` -- a basis of the matroid, assumed to have Python's
      ``frozenset`` interface.
    - ``order`` -- an optional ordering of the base set of M to use
      for the computation.

    Function based on implementation in Matroid._internal.

    TODO:

    Might it be a good idea to make internal and external activity functions
    part of the public interface to the Matroid class?
    """
    N = M.groundset() - B
    A = set()
    order_key = lambda x: x
    if order != None:
        order_key = lambda x: order.index(x) # order should contain all indices
    for e in B:
        if min(M.cocircuit(N | set([e])), key=order_key) == e:
            A.add(e)
    return A


# TODO find a better place to encapsulate this
class ZonUtils(object):
    @staticmethod
    def minimal_sets(subsets):
        """
        Returns the minimal elements (with respect to inclusion) of a collection of subsets.
        """
        min_sets = []
        for s1 in subsets:
            minimal = True
            larger = []
            for s2 in min_sets:
                if set(s1) >= set(s2):
                    minimal = False
                elif set(s1) < set(s2):
                    larger.append(s2)
            for s in larger:
                min_sets.remove(s)
            if minimal:
                min_sets.append(s1)
        return min_sets
