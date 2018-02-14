import sage.matroids.matroid
from sage.combinat.posets.posets import Poset
from sage.combinat.posets.lattices import LatticePoset
from sage.combinat.subset import Subsets

import logging
import sys
logging.basicConfig(stream=sys.stderr, level=logging.WARNING)


class OrderedMatroid(sage.matroids.matroid.Matroid):
    r"""
    Class OrderedMatroid is a class that derives from Matroid.
    It implements functionality induced on a matroid by an ordering of its
    ground set.

    Implements compositing pattern to integrate with built-in Sage Matroid
    functionality.
    """

    def __init__(self, M,  ordered_groundset=None, **kwargs):
        self._ground_matroid = M

        gs = list(M.groundset())
        if ordered_groundset is not None:
            if set(gs) <= set(ordered_groundset):
                gs = sorted(gs, key=(lambda x: ordered_groundset.index(x)))
            else:
                raise ValueError("OrderedMatroid: ordered groundset does not "
                                 "include matroid groundset")
        else:
            gs = sorted(gs, **kwargs)

        def gs_key(x):
            return gs.index(x)

        def gs_cmp(x, y):
            return gs.index(x) - gs.index(y)

        self._gs = gs          # groundset in sorted order
        self._gs_key = gs_key  # groundset key function from this ordered list
        self._gs_cmp = gs_cmp  # groundset cmp function from this ordered list

    def groundset_order(self):
        return list(self._gs)

    def _sorted(self, S, **kwargs):
        return sorted(S, key=self._gs_key, **kwargs)

    def __getattr__(self, item):
        return getattr(self._ground_matroid, item)  # redirection

    def _rank(self, X):
        return self._ground_matroid._rank(X)

    def groundset(self):
        return self._ground_matroid.groundset()

    def _repr_(self):
        S = "Ordered matroid of rank " + str(self.rank()) + " on " \
            + str(self.size()) + " elements"
        return S

    def base_matroid(self):
        """
        Return the underlying unordered matroid of this ordered matroid.
        """
        return self._ground_matroid

    def _circuit_root(self, C):
        return min(C, key=self._gs_key)

    def circuit_root(self, C):
        if self._ground_matroid.is_circuit(C):
            return self._circuit_root(C)
        else:
            raise ValueError(str(C) + " is not a circuit")

    def _cocircuit_root(self, D):
        return min(D, key=self._gs_key)

    def cocircuit_root(self, D):
        if self._ground_matroid.is_cocircuit(D):
            return self._cocircuit_root(D)
        else:
            raise ValueError(str(D) + " is not a cocircuit")

    def _dominant_basis(self, X):
        """
        Return the lex maximal independent subset whose matroid closure
        contains ``X``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The spanning independent subset of ``X`` which is lex maximal with
        respect to the ordering of the matroid.

        EXAMPLES::

            sage: OM = OrderedMatroid(matroids.named_matroids.Vamos())
            sage: sorted( OM._lex_dominant_basis(['c', 'e', 'f', 'g', 'h']) )
            ['c', 'f', 'g', 'h']
            sage: sorted( OM._lex_dominant_basis(['b', 'd', 'e']) )
            ['b', 'd', 'e']
        """
        X_sorted = self._sorted(X, reverse=True)

        I = set()  # independent set
        cl_I = self._ground_matroid.closure(I)
        for x in X_sorted:
            if x not in cl_I:
                I.add(x)
                cl_I = self._ground_matroid.closure(I)

        return frozenset(I)

    def _dominant_cobasis(self, X):
        """
        Return the lex maximal coindependent subset whose matroid coclosure
        contains ``X``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The cospanning coindependent subset of ``X`` which is lex maximal with
        respect to the ordering of the matroid.

        EXAMPLES::

            sage: OM = OrderedMatroid(matroids.named_matroids.Vamos())
            sage: sorted( OM._lex_dominant_cobasis(['c', 'e', 'f', 'g', 'h']) )
            ...
            sage: sorted( OM._lex_dominant_basis(['b', 'd', 'e']) )
            ...
        """
        X_sorted = self._sorted(X, reverse=True)

        cI = set()  # coindependent set
        cl_cI = self._ground_matroid.coclosure(cI)
        for x in X_sorted:
            if x not in cl_cI:
                cI.add(x)
                cl_cI = self._ground_matroid.coclosure(cI)

        return frozenset(cI)

    def _indep_activity(self, I):
        """
        Return the active elements corresponding to an independent set.

        INPUT:

        - ``I`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()`` for which ``self.is_independent(I)``
          returns True.

        OUTPUT:

        A subset of ``self.groundset()`` representing the elements active with
        respect to ``I``.

        EXAMPLES:

            sage: X = Matrix(QQ, [[1, 0], [0, 1], [1, 1], [1, 1]]).transpose()
            sage: M = Matroid(matrix=X)
            sage: OM = OrderedMatroid(M)
            sage: sorted( OM._indep_activity([0, 3]) )
            [2]
            sage: sorted( OM._indep_activity([3]) )
            [2]
            sage: sorted( OM._indep_activity([0,1]) )
            []
            sage: sorted( OM._indep_activity([]) )
            []

            sage: OM = OrderedMatroid(matroids.named_matroids.Vamos())
            sage: sorted( OM._indep_activity(['d', 'e', 'f', 'h']) )
            ['a', 'b', 'c']
        """
        logging.info('called OrderedMatroid._indep_activity for set I = '
                     + str(I))
        logging.debug('sorted groundset self._gs = ' + str(self._gs))

        # output active set
        active = frozenset()

        # independent set ordered from largest to smallest
        indep_elts = sorted(I, key=self._gs_key, reverse=True)

        logging.debug('computing active elements by iterating through flag of'
                      ' flats')

        prev_F = frozenset()
        logging.debug('entering loop')
        for (i, x) in enumerate(indep_elts):
            logging.debug('=== starting loop iteration with (i, x) = '
                          + str((i, x)))

            # flat F spanned by elements x and larger
            F = self._ground_matroid.closure(indep_elts[:i + 1])
            logging.debug('F = ' + str(F))
            logging.debug('prev_F = ' + str(prev_F))

            # restrict to elements which are newly spanned, not in prev_F
            newly_spanned = sorted(F - prev_F, key=self._gs_key)
            logging.debug('newly_spanned = ' + str(newly_spanned))

            # elements smaller than x in this flag component are active wrt I
            new_active = newly_spanned[:newly_spanned.index(x)]
            logging.debug('new_active = ' + str(new_active))
            active |= frozenset(new_active)
            logging.debug('active = ' + str(active))

            # record most recent flat F for next flag computation
            prev_F = F

            logging.debug('=== ending loop iteration')

        return frozenset(active)

    def _coindep_coactivity(self, cI):
        # output coactive set
        coactive = frozenset()

        # independent set ordered from largest to smallest
        coindep_elts = sorted(cI, key=self._gs_key, reverse=True)
        prev_cF = frozenset()
        for (i, x) in enumerate(coindep_elts):
            # coflat cF spanned by elements x and larger
            cF = self._ground_matroid.coclosure(coindep_elts[:i + 1])
            # restrict to elements which are newly spanned, not in prev_cF
            newly_spanned = sorted(cF - prev_cF, key=self._gs_key)
            # elements smaller than x in flag component are coactive wrt cI
            new_coactive = newly_spanned[:newly_spanned.index(x)]
            coactive |= frozenset(new_coactive)
            # record most recent flat F for next flag computation
            prev_cF = cF

        return frozenset(coactive)

    def active_elements(self, X):
        """
        Return the active elements corresponding to an arbitrary subset of the
        matroid ground set.  This is accomplished by finding the lex maximal
        spanning independent set, and computing this independent set's
        corresponding active elements.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A subset of ``self.groundset()`` representing the elements active with
        respect to ``I``.

        EXAMPLES::

            sage: X = Matrix(QQ, [[1, 0], [0, 1], [1, 1], [1, 1]]).transpose()
            sage: M = Matroid(matrix=X)
            sage: OM = OrderedMatroid(M)
            sage: sorted( OM.activity([0, 3]) )
            [2]
            sage: sorted( OM.activity([0, 2, 3]) )
            [2]
            sage: sorted( OM.activity([1, 3]) )
            [0, 2]
            sage: sorted( OM.activity([0, 1, 2, 3]) )
            [0, 2]
        """
        B = self._dominant_basis(X)
        return self._indep_activity(B)

    def coactive_elements(self, X):
        cB = self._dominant_cobasis(X)
        return self._coindep_coactivity(cB)

    def passive_elements(self, X):
        return frozenset(self._gs) - self.active_elements(X)

    def copassive_elements(self, X):
        return frozenset(self._gs) - self.coactive_elements(X)

    def dual(self):
        """
        Returns the dual ordered matroid of ``self``.
        """
        M = self._ground_matroid.dual()
        OM = OrderedMatroid(M, ordered_groundset=self._gs)
        return OM

    def minor(self, contractions=None, deletions=None):
        M = self._ground_matroid.minor(contractions=contractions,
                                       deletions=deletions)
        OM = OrderedMatroid(M, ordered_groundset=self._gs)
        return OM

    def delete(self, X):
        return self.minor(deletions=X)

    def contract(self, X):
        return self.minor(contractions=X)

    def _set_to_str(self, X):
        elts = self._sorted(X)
        elts = [str(x) for x in elts]
        return ''.join(elts)

    def external_order(self):
        M = self._ground_matroid
        passives = {I: self.passive_elements(I) - I
                    for I in M.independent_sets()}

        def label(I):
            # str_I = self._set_to_str(I)
            str_EP = self._set_to_str(passives[I])
            return str_EP  # + " (" + str_I + ")"
        labels = {passives[I]: label(I) for I in passives}

        def ext_cmp(S, T):
            return S.issubset(T)

        P = Poset(data=(passives.values(), ext_cmp), element_labels=labels)
        return P

    def internal_order(self):
        return self.dual().external_order()

    def internal_external_order(self):
        E = frozenset(self._gs)

        passives = {}
        for S in Subsets(E):
            S = frozenset(S)
            EP = self.passive_elements(S) - S
            IP = self.copassive_elements(E - S) & S
            passives[S] = (EP, IP)

        # TODO poset fails to be a lattice in simple cases... Is the
        # construction wrong?

        # for S in Subsets(E):
        #     S = frozenset(S)
        #     EP = self.active_elements(S)
        #     IP = self.coactive_elements(E - S)
        #     passives[S] = (EP, IP)
        #     print self._set_to_str(S)
        #     print "a" + self._set_to_str(EP), "b" + self._set_to_str(IP)

        labels = {}
        for S in passives:
            pS = passives[S]
            str_ext = self._set_to_str(passives[S][0])
            str_int = self._set_to_str(passives[S][1])
            labels[pS] = str_ext + "/" + str_int

        def int_ext_cmp(S, T):
            return S[0].issubset(T[0]) and T[1].issubset(S[1])

        P = LatticePoset(data=(passives.values(), int_ext_cmp),
                         element_labels=labels)
        return P
