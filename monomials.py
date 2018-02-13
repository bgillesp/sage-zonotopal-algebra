from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.functions.other import binomial
from sage.rings.polynomial.polydict import ETuple
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.structure.parent import Parent
from sage.structure.parent import Set_generic
from sage.structure.unique_representation import UniqueRepresentation


class Monomials(Set_generic, UniqueRepresentation):
    r"""
    Set representing a collection of monomials coming from a polynomail ring.

    INPUT:

    - ``P`` -- an ambient polynomial ring
    - ``degree`` (default: ``(0, +Infinity)``) an integer representing the
      single fixed degree of monomials, or a tuple of length 2 consisting of
      lower and upper bounds (inclusive) for the degree of monomials

    OUTPUT:

    - return a Monomials object with given ambient polynomial ring and degrees

    EXAMPLES:

        sage: P.<x, y> = PolynomialRing(QQ)
        sage: monoms = Monomials(P, (0, 2))
        sage: [m for m in monoms]
        [1, x, y, x^2, x*y, y^2]
        sage: x^2 in monoms
        True
        sage: x + y in monoms
        False
    """

    def __init__(self, P, degree=(0, +Infinity)):
        r"""
        The set of monomials of a polynomial ring.
        P is the polynomial ring.
        degree is either a tuple with the minimum and maximum degrees,
        inclusive, or an integer with the single fixed degree
        """
        # is there a good way to check that P is a polynomial ring?
        self._poly_ring = P
        if (type(degree) != tuple) and (degree not in ZZ):
            raise ValueError("degree must be a nonnegative integer or a tuple")
        if degree in ZZ:
            self._min_deg = degree
            self._max_deg = degree
        else:
            self._min_deg = degree[0]
            self._max_deg = degree[1]
        m, M = self._min_deg, self._max_deg
        if (m not in NN) or (M not in NN and M != +Infinity):
            raise ValueError("degree bounds must be nonnegative integers,"
                             " or +Infinity for upper bound")
        if self._max_deg < +Infinity:
            super(Monomials, self).__init__(
                facade=self._poly_ring, category=EnumeratedSets.Finite())
        else:
            super(Monomials, self).__init__(
                facade=self._poly_ring, category=EnumeratedSets.Infinite())

    def min_deg(self):
        return self._min_deg

    def max_deg(self):
        return self._max_deg

    def _deg_string(self):
        m, M = self._min_deg, self._max_deg
        if M == m:
            return " with degree %d" % m
        elif m > 0 and M < +Infinity:
            return " with degree at least %d and at most %d" % (m, M)
        elif m > 0:
            return " with degree at least %d" % m
        elif M < +Infinity:
            return " with degree at most %d" % M
        else:
            return ""

    def _repr_(self):
        return "Monomials%s in %s" % (self._deg_string(), self._poly_ring)

    def __contains__(self, elt):
        p = self._poly_ring(elt)
        deg = p.degree()
        return (p.is_monomial()
                and deg >= self._min_deg
                and deg <= self._max_deg)

    def _lex_iterator(self, d):
        length = self._poly_ring.ngens()
        return IntegerListsLex(n=d, length=length, element_constructor=list)

    def __iter__(self):
        deg = self._min_deg
        # iterate through individual degrees, small to large
        while deg <= self._max_deg:
            it = iter(self._lex_iterator(deg))
            try:
                while True:
                    yield self._monomial_from_degrees(it.next())
            except GeneratorExit as e:
                raise e
            # handle termination of intermediate inverse lex iterators
            except StopIteration:
                pass
            deg += 1

    def an_element(self):
        it = iter(self)
        return it.next()

    def _monomial_from_degrees(self, v):
        P = self._poly_ring
        gens = P.gens()
        p = P.one()
        for i, t in enumerate(v):
            p *= gens[i]**t
        return p

    # TODO potentially implement a next(self, elt) method to compute more
    # efficiently than the default EnumeratedSets implementation

    def cardinality(self):
        if self._max_deg == Infinity:
            return Infinity
        else:
            n = self._poly_ring.ngens()
            k1, k2 = self._min_deg, self._max_deg
            return binomial(n + k2, n) - binomial(n + k1 - 1, n)

    def _element_constructor_(self, elt):
        r"""
        Construct a monomial in the appropriate polynomial ring, if it lies
        in the correct range of degrees.
        """
        if elt in self:
            return self._poly_ring(elt)
        else:
            raise ValueError(
                "Value %s is not a monomial%s" % (elt, self._deg_string()))

    def __getitem__(self, index):
        # if a tuple of exponents
        if isinstance(index, tuple) or isinstance(index, ETuple):
            mon = self._monomial_from_degrees(index)
            return self._element_constructor_(mon)
        else:
            return Parent.__getitem__(self, index)
