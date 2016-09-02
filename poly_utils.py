from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.free_module import CombinatorialFreeModuleElement
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.functions.other import factorial
from sage.matrix.constructor import Matrix
from sage.misc.misc_c import prod
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.sets.set import Set
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.parent import Set_generic
from sage.structure.unique_representation import UniqueRepresentation

class Monomials(Set_generic, UniqueRepresentation):
    def __init__(self, P, degree=(0, +Infinity)):
        r"""
        The set of monomials of a polynomial ring.
        P is the polynomial ring
        degree is either a tuple with the minimum and maximum degrees, inclusive,
        or an integer with the single fixed degree
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
            raise ValueError("degree bounds must be nonnegative integers, or +Infinity for upper bound")
        if self._max_deg < +Infinity:
            super(Monomials, self).__init__(facade = self._poly_ring, category = EnumeratedSets.Finite())
        else:
            super(Monomials, self).__init__(facade = self._poly_ring, category = EnumeratedSets.Infinite())

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
        m, M = self._min_deg, self._max_deg
        return "Monomials%s in %s" % (self._deg_string(), self._poly_ring)

    def __contains__(self, elt):
        p = self._poly_ring(elt)
        deg = p.degree()
        return p.is_monomial() and deg >= self._min_deg and deg <= self._max_deg

    def _lex_iterator(self, d):
        m, M = self._min_deg, self._max_deg
        l = self._poly_ring.ngens()
        return IntegerListsLex(n = d, length = l, element_constructor = list)

    def __iter__(self):
        P = self._poly_ring
        deg = self._min_deg
        # iterate through individual degrees, small to large
        while deg <= self._max_deg:
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
            return binomial(n + k2, n) - binomial(n + k1 - 1, n)

    def _element_constructor_(self, elt):
        r"""
        Construct a monomial in the appropriate polynomial ring, if it lies
        in the correct range of degrees.
        """
        if elt in self:
            return self._poly_ring(elt)
        else:
            raise ValueError("Value %s is not a monomial%s" % (elt, self._deg_string()))

    def __getitem__(self, index):
        if type(index) is tuple:
            mon = self._poly_ring.monomial( *index )
            return self._element_constructor_(mon)
        else:
            return Parent.__getitem__(self, index)

# very small extension of CombinatorialFreeModuleElement which adds the ability
# to convert an element back to the appropriate polynomial class
class PolynomialFreeModuleElement(CombinatorialFreeModuleElement):
    def to_polynomial(self):
        terms = [coeff*p for p, coeff in self.monomial_coefficients().items()]
        return sum(terms)

class PolynomialFreeModule(CombinatorialFreeModule, UniqueRepresentation, Parent):
    def __init__(self, P, basis=None):
        r"""
        Warning: for infinite user specified bases, doesn't check that the bases
        are linearly independent
        """
        self._poly_ring = P
        if basis is None:
            basis = Monomials(self._poly_ring)
        if type(basis) == type(Monomials(self._poly_ring)): # infinite monomials basis
            self._basis = basis
            self._dimn = Infinity
            self._name = 'Module generated by all monomials%s in %s' \
                % (self._basis._deg_string(), self._poly_ring)
            self._converter = self._MonomialsConverter(self._poly_ring, self, self._basis)
        elif type(basis) == type(Monomials(self._poly_ring, (0, 1))): # finite monomials basis
            # cast basis to a list here to play nice
            self._basis = basis
            self._dimn = len(self._basis)
            self._name = 'Module generated by all monomials%s in %s' \
                % (self._basis._deg_string(), self._poly_ring)
            # in the converter, use the Monomials class rather than a list for efficiency
            self._converter = self._MonomialsConverter(self._poly_ring, self, self._basis)
        elif Set(basis).cardinality() < Infinity: # other finite basis
            # cast basis to a list here because lazy loading doesn't play nice
            # with some of methods for CombinatorialFreeModule, especially gens()
            self._basis = list(basis)
            self._dimn = len(self._basis)
            if all( [self._poly_ring(p).is_monomial() for p in self._basis] ): # all are monomials
                self._name = 'Module generated by basis of %d monomials in %s' \
                    % (self._dimn, self._poly_ring)
                self._converter = self._MonomialsConverter(self._poly_ring, self, self._basis)
            else: # general finite list of polynomials
                self._name = 'Module generated by basis of %d polynomials in %s' \
                    % (self._dimn, self._poly_ring)
                self._converter = self._FiniteBasisConverter(self._poly_ring, self, self._basis)
        else:  # general infinite basis
            self._basis = basis
            self._dimn = Infinity
            self._name = 'Module generated by infinite basis of polynomials in %s' \
                % self._poly_ring
            self._converter = self._InfiniteBasisConverter(self._poly_ring, self, self._basis)

        super(PolynomialFreeModule, self).__init__(R=self._poly_ring.base_ring(), basis_keys=self._basis)

        # Set up formatting
        bracket = ['(', ')']
        self.print_options(prefix="", bracket=bracket)

    # TODO what the heck does self.basis() return??

    class _MonomialsConverter:
        def __init__(self, P, comb_mod, basis):
            self._poly_ring = P
            self._module = comb_mod
            self._basis = basis

            # assert: all elements in basis are monomials
            if type(self._basis) == list:
                if len(self._basis) != len(set(self._basis)):
                    raise ValueError("Basis polynomials are not linearly independent")

        def convert(self, p):
            module_p = self._module.zero()
            for exponent, coeff in p.dict().items():
                monom = self._poly_ring.monomial(*exponent)
                if not monom in self._basis:
                    raise ValueError("Value %s is not spanned by the basis polynomials" % p)
                module_p += coeff*self._module.monomial( monom )
            return module_p


    class _FiniteBasisConverter:
        def __init__(self, P, comb_mod, basis):
            r"""
            basis should be a finite set of polynomials
            """
            self._poly_ring = P
            self._module = comb_mod
            self._basis = basis

            max_deg = max([self._poly_ring(b).degree() for b in self._basis])
            monoms = []
            for b in self._basis:
                b = self._poly_ring(b)
                monoms += b.monomials()
            monoms_list = tuple(Set(monoms))

            # check if the basis represented in terms of Monomials is efficient
            degs = [self._poly_ring(m).degree() for m in monoms]
            min_deg, max_deg = min(degs), max(degs)
            monoms_obj = Monomials(self._poly_ring, (min_deg, max_deg+1))

            if monoms_obj.cardinality() < 2*len(monoms_list):
                computational_basis = monoms_obj
            else:
                computational_basis = monoms_list
            self._monomial_module = PolynomialModule(P=self._poly_ring, basis=computational_basis)
            self._basis_mat = \
                Matrix([self._monomial_module(b).to_vector() for b in self._basis]).transpose()
            if self._basis_mat.ncols() > self._basis_mat.rank():
                raise ValueError("Basis polynomials are not linearly independent")

        def convert(self, p):
            try:
                p_vect = self._monomial_module(p).to_vector()
                decomp = self._basis_mat.solve_right(p_vect)
            except ValueError:
                raise ValueError("Value %s is not spanned by the basis polynomials" % p)
            polys = [v[1]*self._module.monomial( v[0] ) for v in zip(self._basis, decomp)]
            module_p = sum(polys, self._module.zero())
            return module_p

    class _InfiniteBasisConverter:
        def __init__(self, P, comb_mod, basis):
            self._poly_ring = P
            self._module = comb_mod
            self._basis = basis
            # TODO implement some sort of caching for intermediate objects for conversion

        def convert(self, p):
            # TODO check if polynomial is contained in span of basis in general
            # for now assume that infinite basis is in increasing order of degree,
            # and collect basis elements until the degree exceeds that of input poly
            # Note: This is not necessarily a good heuristic without further processing
            # of the basis beforehand; always works for homogeneous polynomials though
            # TODO for now, what characteristics of the basis guarantee correctness here?
            deg = p.degree()
            fin_basis = []
            it = iter(self._basis)
            b = it.next()
            while b.degree() <= deg:
                fin_basis.append(b)
                b = it.next()

            C = PolynomialModule._FiniteBasisConverter(self._poly_ring, self._module, fin_basis)
            return C.convert(p)

    def _element_constructor_(self, elt):
        return self.convert(elt)

    def convert(self, p):
        if not p in self._poly_ring:
            raise ValueError("Value %s is not a polynomial in %s" % (p, self._poly_ring))
        # ensure that the type of p is actually a polynomial, not, e.g., an integer
        p = self._poly_ring(p)
        return self._converter.convert(p)

    Element = PolynomialFreeModuleElement

# TODO allow conversions between elements of different modules over the same
# polynomial space, and implement conversions/coercions

# would like to be able to:
# * convert between spaces seamlessly, find join and meet of spaces,
# * represent a basis of a polynomial space (as in an ideal with generators)
# * transition easily between representations in different bases

# perhaps it can automatically transition dimension, e.g. cast things up in dimension
# automatically but trim dimension only by explicit request


class PolyUtils:
    @staticmethod
    def poly_deriv(p,q):
        g = p.parent().gens();
        s = 0;
        for e_tup, coeff in p.dict().iteritems():
            diff_list = [];
            for v,e in zip(g,e_tup):
                diff_list.extend([v]*e);
            s += coeff*q.derivative(diff_list);
        return s

    @staticmethod
    def diff_bilinear_form(p,q):
        """
        The bilinear form on polynomials derived from the inner product on
        forms, corresponding to the differential operation ``p(D)q`` evaluated
        at ``x = 0``.
        """
        n_vars = len(p.parent().gens())
        zero = [0]*n_vars
        return (PolyUtils.poly_deriv(p,q))(zero)

    # TODO do we need this method given that we have the Monomials class?
    @staticmethod
    def monomial_from_degrees(P,v):
        gens = P.gens()
        p = P.one()
        for i,t in enumerate(v):
            p *= gens[i]**t
        return p

    @staticmethod
    def linear_form(P, vec):
        terms = [coeff*var for coeff, var in zip(vec, P.gens())]
        return sum(terms, P.zero())

    @staticmethod
    def pure_tensor(P, vects, indices=None):
        """
        Form a product of linear forms using vectors in the list `vects` for the
        forms and the elements of `indices` as the indices.
        Items in `indices` may be repeated
        """
        if indices == None:
            indices = range(len(vects))
        terms = [PolyUtils.linear_form(P, vects[i]) for i in indices]
        return prod(terms, P.one())

    @staticmethod
    def poly_dual_basis(P, poly_basis):
        # compute a dual basis for an input *homogeneous* basis
        deg = max([p.degree() for p in poly_basis])
        poly_module = PolynomialFreeModule( P, basis=Monomials(P, (0, deg+1)) )
        bilinear_form_coeffs = []
        for b in poly_module.basis().keys():
            # each b is a monomial in self.Pi of degree at most deg
            b = P(b)
            bilinear_form_coeffs.append( prod( map(factorial, b.degrees()) ) )
        A = Matrix([poly_module(p).to_vector() for p in poly_basis])
        D = Matrix.diagonal(bilinear_form_coeffs, sparse=False)
        B = (A*D*A.transpose()).inverse()
        dual_basis = []
        for col in B.columns():
            q = sum( [coeff*p for coeff,p in zip(col, poly_basis)] )
            dual_basis.append(q)
        return dual_basis
