from set_utils import SetUtils
from sage.misc.misc import powerset
from sage.combinat.posets.posets import Poset


# TODO M.is_coindependent is probably more efficient than iterating through
# bases individually as in the below general implementations.  More efficient
# methods in general?
# def short_subsets(self):
#     M = self._matroid()
#     for s in powerset(M.groundset()):
#         if M.is_coindependent(s):
#             yield s
#
# def long_subsets(self):
#     M = self._matroid()
#     for s in powerset(M.groundset()):
#         if not M.is_coindependent(s):
#             yield s

# TODO implement more efficient method for determining generalized cocircuits,
# which are the minimal long subsets
def long_subsets(M, int_subsets=None):
    r"""
    Generates the long subsets of a matroid ``M`` with respect to a collection
    ``int_subsets`` of intersecting subsets. These are the subsets of ``M``
    which have nontrivial intersection with every set in ``int_subsets``.  This
    is the complement of the set of short subsets.

    INPUT:

    - ``M`` -- matroid object for which to compute long subsets

    - ``int_subsets`` -- a list of bases (default: None) in ``M`` to be used for
      the computation.  If None, then uses the collection of bases of ``M``, the
      classical definition for long subsets.

    OUTPUT:

    A generator for the long subsets of ``M`` with respect to the given
    intersecting subsets.

    EXAMPLES:

    A simple matroid gives the following long sets with respect to the matroid's
    bases:

    ::

        sage: X = Matrix(QQ, [[1, 0], [0, 1], [1, 1]]).transpose(); X
        [1 0 1]
        [0 1 1]
        sage: M = Matroid(X)
        sage: for l in long_subsets(M):
        ....:     print l
        ....:
        [0, 1]
        [0, 2]
        [1, 2]
        [0, 1, 2]

    Further specifying a particular collection of intersecting subsets:

    ::

        int_subsets = [[0, 1], [2]]
        sage: for l in long_subsets(M, int_subsets):
        ....:     print l
        ....:
        [0, 2]
        [1, 2]
        [0, 1, 2]
    """
    if int_subsets == None:
        int_subsets = M.bases()
    int_subsets = list(int_subsets)
    for s in powerset(M.groundset()):
        for b in int_subsets:
            if len( set(b).intersection(set(s)) ) == 0:
                break
        else:
            yield s

def short_subsets(M, int_subsets=None):
    r"""
    Generates the short subsets of a matroid ``M`` with respect to a collection
    ``int_subsets`` of intersecting subsets. These are the subsets of ``M``
    which avoid some set in ``int_subsets``.  This is the complement of the set
    of long subsets.

    INPUT:

    - ``M`` -- matroid object for which to compute long subsets

    - ``int_subsets`` -- a list of bases (default: None) in ``M`` to be used for
      the computation.  If None, then uses the collection of bases of ``M``, the
      classical definition for long subsets.

    OUTPUT:

    A generator for the long subsets of ``M`` with respect to the given
    intersecting subsets.

    EXAMPLES:

    A simple matroid gives the following short sets with respect to the
    matroid's bases:

    ::

        sage: X = Matrix(QQ, [[1, 0], [0, 1], [1, 1]]).transpose(); X
        [1 0 1]
        [0 1 1]
        sage: M = Matroid(X)
        sage: for l in short_subsets(M):
        ....:     print l
        ....:
        []
        [0]
        [1]
        [2]

    Further specifying a particular collection of intersecting subsets:

    ::

        int_subsets = [[0, 1], [2]]
        sage: for l in short_subsets(M, int_subsets):
        ....:     print l
        ....:
        []
        [0]
        [1]
        [0, 1]
        [2]
    """
    if int_subsets == None:
        int_subsets = M.bases()
    int_subsets = list(int_subsets)
    for s in powerset(M.groundset()):
        for b in int_subsets:
            if len( set(b).intersection(set(s)) ) == 0:
                yield s
                break

def cocircuits(M, bases=None):
    if bases == None:
        bases = M.bases()
    long_sets = long_subsets(M, bases)
    return SetUtils.minimal_sets(long_sets)

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

def external_order_leq(M, order=None):
    def leq(B1, B2):
        B1 = frozenset(B1)
        B2 = frozenset(B2)
        if not (M.is_basis(B1) and M.is_basis(B2)):
            raise ValueError("Inputs are not both bases")
        # B1 is ext leq B2 if B2 is contained in B1 union its ext active elts
        return B2.issubset( B1.union( matroid_external(M, B1, order) ) )
    return leq

def external_order_poset(M, order=None):
    leq = external_order_leq(M, order=order)
    bases = [tuple(b) for b in M.bases()]
    P = Poset((bases, leq))
    return P
