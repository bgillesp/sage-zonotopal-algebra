import poly_utils
from set_utils import SetUtils
from matroid_utils import external_order_leq
from sage.misc.cachefunc import cached_method
from abstract_zonotopal_algebra import AbstractZonotopalAlgebra




class ForwardExchangeZonotopalAlgebra(AbstractZonotopalAlgebra):
    def __init__(self,X,bases,varNames="x"):
        AbstractZonotopalAlgebra.__init__(self, X, varNames)

        # initialize forward exchange bases by finding maximal elements
        M = self._matroid()
        # our external order is based on the external activity for the reverse ordering
        reverse_order = range(X.ncols() - 1, -1, -1)
        self.ext_order = external_order_leq(M, order=reverse_order)
        # find the maximal bases
        self.max_bases = SetUtils.max_elements(bases, self.ext_order)


    # TODO potentially implement a poset class corresponding to the external order

    def __repr__(self):
        return "Forward Exchange Zonotopal Algebra over " + str(self.base_field()) + " with matrix\n" + str(self.matrix())

    # TODO static
    def _cocircuit_basis_order(self, cocirc, basis):
        boundary = min( set(cocirc).intersection(set(basis)) )
        truncated_elts = filter((lambda x: x <= boundary), cocirc)
        return len(truncated_elts)

    def _cocircuit_order(self, cocirc):
        orders = [self._cocircuit_basis_order(cocirc, B) for B in self.max_bases]
        return max(orders)

    # TODO static
    def _truncate_ordered(self, set, size):
        l = list(set)
        l.sort()
        return l[:size]

    def _truncated_cocircuit(self, cocirc):
        order = self._cocircuit_order(cocirc)
        return self._truncate_ordered(cocirc, order)

    def _truncated_cocircuits(self):
        for cocirc in self._cocircuits():
            yield self._truncated_cocircuit(cocirc)

    @cached_method
    def I_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        n = self._matroid().size()
        for normal,cocirc in zip( self._facet_hyperplane_normals(), self._cocircuits() ):
            # degree of linear form is the size of the truncated cocircuit
            deg = self._cocircuit_order(cocirc)
            gen = poly_utils.linear_form(P, normal)**deg
            gens.append(gen)
        return gens

    @cached_method
    def J_ideal_gens(self):
        gens = []
        P = self.polynomial_ring()
        X_cols = self.matrix().columns()
        for cocirc in self._cocircuits():
            fe_cocirc = self._truncated_cocircuit(cocirc)
            gen = poly_utils.pure_tensor(P, X_cols, fe_cocirc)
            gens.append(gen)
        return gens
