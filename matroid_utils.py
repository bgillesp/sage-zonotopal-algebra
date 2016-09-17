from set_utils import SetUtils
from sage.misc.misc import powerset


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
def long_subsets(M, bases=None):
    if bases == None:
        bases = M.bases()
    bases = list(bases)
    for s in powerset(M.groundset()):
        for b in bases:
            if len( set(b).intersection(set(s)) ) == 0:
                break
        else:
            yield s

def short_subsets(M, bases=None):
    if bases == None:
        bases = M.bases()
    bases = list(bases)
    for s in powerset(M.groundset()):
        s = Set(s)
        for b in bases:
            b = Set(b)
            if len( b.intersection(s) ) == 0:
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