
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.categories.enumerated_sets import EnumeratedSets
#from sage.rings.integer import Integer



class Monomials(UniqueRepresentation, Parent):
    def __init__(self, P, degs=(0, +Infinity)):
        r"""
        The set of monomials of a polynomial ring.
        P is the polynomial ring
        degs is the minimum degree inclusive and the maximum degree exclusive
        """
        # is there a good way to check that P is a polynomial ring?
        self._poly_ring = P
        #self.Element = MonomialElement
        self._min_deg = degs[0]
        self._max_deg = degs[1]
        if self._max_deg < +Infinity:
            Parent.__init__(self, category = EnumeratedSets.Finite())
        else:
            Parent.__init__(self, category = EnumeratedSets.Infinite())

        self._populate_coercion_lists_( embedding=MonomialsMorphism(self, self._poly_ring) )

    def _deg_string(self):
        m, M = self._min_deg, self._max_deg - 1

        if m > 0 and M < +Infinity:
            return " with degree at least %d and at most %d" % (m, M)
        elif m > 0:
            return " with degree at least %d" % m
        elif M < +Infinity:
            return " with degree at most %d" % M
        else:
            return ""

    def _repr_(self):
        m, M = self._min_deg, self._max_deg - 1
        return "The enumerated set of monomials%s in %s" % \
            (self._deg_string(), P)

        # if m > 0 and M < +Infinity:
        #     return "%s with degree at least %d and at most %d" \
        #         % (base_str, self._min_deg, self._max_deg - 1)
        # elif m > 0:
        #     return "%s with degree at least %d" % (base_str, self._min_deg)
        # elif m < +Infinity:
        #     return "%s with degree at most %d" % (base_str, self._max_deg)
        # else:
        #     return base_str

    def __contains__(self, elt):
        p = self._poly_ring(elt)
        deg = p.degree()
        return p.is_monomial() and deg >= self._min_deg and deg < self._max_deg

    def _lex_iterator(self, d):
        m, M = self._min_deg, self._max_deg - 1
        l = self._poly_ring.ngens()
        return IntegerListsLex(n = d, length = l, element_constructor = list)

    def __iter__(self):
        P = self._poly_ring
        deg = self._min_deg
        # iterate through individual degrees, small to large
        while deg < self._max_deg:
            it = iter( self._lex_iterator(deg) )
            continueBlock = False
            try:
                while True:
                    yield P.monomial( *it.next() )
            except GeneratorExit as e:
                raise e
            except: # handle when intermediate inverse lexicographic iterators end
                pass
            deg += 1

    def an_element(self):
        it = iter(self)
        return it.next()

    # TODO potentially implement a next(self, elt) method to compute more
    #  more efficiently than the default EnumeratedSets implementation

    def cardinality(self):
        if self._max_deg == Infinity:
            return Infinity
        else:
            n = self._poly_ring.ngens()
            k1, k2 = self._min_deg, self._max_deg
            return binomial(n + k2 - 1, n) - binomial(n + k1 - 1, n)

    # def _element_constructor_(self, elt):
    #     r"""
    #     Construct a monomial in the appropriate polynomial ring, if it lies
    #     in the correct range of degrees.
    #     """
    #     if elt in self:
    #         return self._poly_ring(elt)
    #     else:
    #         raise ValueError("Value %s is not a monomial%s" % (elt, self._deg_string()))

    def __getitem__(self, index):
        if type(index) is tuple:
            mon = self._poly_ring.monomial( *index )
            return self._element_constructor_(mon)
        else:
            return Parent.__getitem__(self, index)

    Element = MonomialElement

class MonomialElement(RingElement):
    def __init__(self, parent, p):
        Element.__init__(self, parent)
        if p in parent:
            self._p = parent._poly_ring(p)
        else:
            raise ValueError("Value %s is not a monomial%s" % (p, parent._deg_string()))

    def _repr_(self):
        return self._p._repr_()

    def _value(self):
        return self._p

class MonomialsMorphism(Morphism):
    def __init__(self, monom_set, poly_ring):
        Morphism.__init__(self, monom_set, poly_ring)

    def _call_(self, x):
        return self.domain()._poly_ring(x._value())

class PolynomialModule(CombinatorialFreeModule, UniqueRepresentation, Parent):
    def __init__(self, P, basis=None):
        r"""
        Warning: for infinite user specified bases, doesn't check that the bases
        are linearly independent
        """
        self._poly_ring = P
        self._basis = basis
        if not basis:
            basis = Monomials(P)
            self._dimn = Infinity
            self._uses_brackets = True
            self._name = 'Module of polynomials in %s with monomial basis' % P
        else:
            self._dimn = Set(self._basis).cardinality()
            self._uses_brackets = True
            if self._dimn == Infinity():
                self._name = 'Module of polynomials in %s with infinite basis' % P
            else:
                # TODO check that the given basis is linearly independent
                self._name = 'Module of polynomials in %s with basis of size %d' % (P, self._dimn)

        CombinatorialFreeModule.__init__(self, R=self._poly_ring.base_ring(),\
            basis_keys=self._basis)

        # Set up formatting
        bracket = ['(', ')'] if self._uses_brackets else False
        self.print_options(prefix="m", bracket=bracket)

class FiniteDegreePolynomialModule(CombinatorialFreeModule, UniqueRepresentation, Parent):
    def __init__(self, P, max_deg=0):
        self._poly_ring = P
        self._max_deg = max_deg
        self._basis = Monomials(self._poly_ring, (0, self._max_deg + 1))
        self._name = 'Module of polynomials of degree at most %d in %s with monomial basis' \
            % (self._max_deg, self._poly_ring)
        self.print_options(prefix="m", bracket=False)

        CombinatorialFreeModule.__init__(self, self._poly_ring.base_ring(), self._basis)

    def __call__(self, o):
        return self._element_constructor_(o)

    def _convert_polynomial_(self, p):
        poly_items = [(poly_vectors(exponent), coeff) for exponent, coeff in p.dict().items()]

    def _element_constructor_(self, elt):
        if not elt in self._poly_ring or elt.degree() > self._max_deg:
            raise ValueError("Value %s is not a polynomial of degree at most %d in %s" \
                % (elt, self._max_deg, self._poly_ring))
        else:
            b = self.basis()
            poly_items = [(poly_vectors(exponent), coeff) for exponent, coeff in p.dict().items()]


class OldPolynomialModule(CombinatorialFreeModule):
    def __init__(self, P, k):
        self._name = 'Module of degree at most %d polynomials' % k
        CombinatorialFreeModule.__init__(self, P.base_ring(), Monomials(P, (0,k+1)))
        self.print_options(prefix="m", bracket=["(",")"])

#

# IntegerListsLex(max_sum=n, length=k)

# need to be able to convert between spaces seamlessly, find join and meet of spaces,
# decompose a polynomial to a vector form in terms of an arbitrary basis
# reconstitute a polynomial from a vector in terms of an arbitrary basis

# convert elements to polynomials

# give canonical polynomial space for this vector space

# convert polynomials to element of free module

# represent a basis of a polynomial space (as in an ideal with generators)

# transition between representations in different bases

# polynomial vector space represents a vector space of polynomials, needs to be able to
# represent arbitrary polynomials in an appropriate finite dimensional vector space

# perhaps it can automatically transition dimension, e.g. cast things up in dimension
# automatically but trim dimension only by explicit request

# define a class of finite polynomial vector space
