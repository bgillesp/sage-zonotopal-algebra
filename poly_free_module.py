# TODO remove the following (included for testing purposes outside sage):
import sys
sys.path.insert(0, '/home/bgillespie/Desktop/software-dev/zonotopal-algebra/'
                   'sage-zonotopal-algebra/lib')

from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.constructor import Matrix
from sage.rings.infinity import Infinity
from sage.sets.set import Set
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from monomials import Monomials


class PolynomialFreeModuleElement(CombinatorialFreeModule.Element):
    r"""
    Element class for PolynomialFreeModule.

    The class is a very small extension of CombinatorialFreeModuleElement
    which adds the ability to convert an element back to the appropriate
    polynomial class.
    """
    def to_polynomial(self):
        monom_coeffs = self.monomial_coefficients()
        terms = [coeff * p for p, coeff in monom_coeffs.items()]
        return sum(terms)


class PolynomialFreeModule(CombinatorialFreeModule, UniqueRepresentation,
                           Parent):
    r"""
    Class PolynomialFreeModule implements a combinatorial free module
    whose underlying elements come from a polynomial ring.  This structure
    is useful when considering the induced module or vector space
    structure of a polynomial ring coming from the multiplication
    of the coefficient field or ring.

    INPUT:

    - ``P`` -- an ambient polynomial ring
    - ``basis`` -- (default: ``None``) a collection of polynomials in P
      to use as a basis for the resulting module.  This must be a hashable type
      in order to work properly with the ``UniqueRepresentation`` class, e.g. a
      ``frozenset`` object or a ``tuple``.

    OUTPUT:

    - return a PolynomialFreeModule

    EXAMPLES:

        sage: P.<x> = PolynomialRing(QQ)
        sage: basis = (1, x, x^2 + x)
        sage: M = PolynomialFreeModule(P, basis)
        sage: M(1 + x^2)
        (1) - (x) + (x^2 + x)
        sage: M(1 + x^2).to_vector()
        (1, -1, 1)
        sage: M(1 + x^3)
        Traceback (most recent call last):
        ...
        ValueError: Value x^3 + 1 is not spanned by the basis polynomials
    """

    def __init__(self, P, basis=None):
        r"""
        Warning: for infinite user specified bases, doesn't check that the
        bases are linearly independent
        """
        self._poly_ring = P
        if basis is None:
            basis = Monomials(self._poly_ring)
        # infinite monomials basis
        if isinstance(basis, type(Monomials(self._poly_ring))):
            self._basis = basis
            self._dimn = Infinity
            self._name = 'Module generated by all monomials%s in %s' \
                % (self._basis._deg_string(), self._poly_ring)
            self._converter = self._MonomialsConverter(
                self._poly_ring, self, self._basis)
        # finite monomials basis
        elif isinstance(basis, type(Monomials(self._poly_ring, (0, 1)))):
            # cast basis to a list here to play nice
            self._basis = basis
            self._dimn = len(self._basis)
            self._name = (
                'Module generated by all monomials%s in %s'
                % (self._basis._deg_string(), self._poly_ring))
            # use the Monomials class rather than a list for efficiency:
            self._converter = self._MonomialsConverter(
                self._poly_ring, self, self._basis)
        # other finite basis
        elif Set(basis).cardinality() < Infinity:
            # cast basis to a list because lazy loading doesn't play nice with
            # some of methods for CombinatorialFreeModule, especially gens()
            self._basis = tuple(basis)
            self._dimn = len(self._basis)
            # all are monomials:
            if all([self._poly_ring(p).is_monomial() for p in self._basis]):
                self._name = (
                    'Module generated by basis of %d monomials in %s'
                    % (self._dimn, self._poly_ring))
                self._converter = self._MonomialsConverter(
                    self._poly_ring, self, self._basis)
            # general finite list of polynomials
            else:
                self._name = (
                    'Module generated by basis of %d polynomials in %s'
                    % (self._dimn, self._poly_ring))
                self._converter = self._FiniteBasisConverter(
                    self._poly_ring, self, self._basis)
        else:  # general infinite basis
            self._basis = basis
            self._dimn = Infinity
            self._name = (
                'Module generated by infinite basis of polynomials in %s'
                % self._poly_ring)
            self._converter = self._InfiniteBasisConverter(
                self._poly_ring, self, self._basis)

        super(PolynomialFreeModule, self).__init__(
            R=self._poly_ring.base_ring(), basis_keys=self._basis)

        # Set up formatting
        bracket = ['(', ')']
        self.print_options(prefix="", bracket=bracket)

    # TODO check what self.basis() returns for CombinatorialFreeModule objects

    class _MonomialsConverter:
        def __init__(self, P, comb_mod, basis):
            self._poly_ring = P
            self._monoms = Monomials(P)
            self._module = comb_mod
            self._basis = basis

            # assert: all elements in basis are monomials
            if type(self._basis) == list:
                if len(self._basis) != len(set(self._basis)):
                    raise ValueError(
                        "Basis polynomials are not linearly independent")

        def convert(self, p):
            module_p = self._module.zero()
            for exponent, coeff in p.dict().items():
                monom = self._monoms[exponent]
                if monom not in self._basis:
                    raise ValueError(
                        "Value %s is not spanned by the basis polynomials" % p)
                module_p += coeff * self._module.monomial(monom)
            return module_p

    class _FiniteBasisConverter:
        def __init__(self, P, comb_mod, basis):
            r"""
            Basis should be a finite set of polynomials
            """
            self._poly_ring = P
            self._module = comb_mod
            self._basis = basis

            max_deg = max([self._poly_ring(b).degree() for b in self._basis])
            monoms = []
            for b in self._basis:
                poly = self._poly_ring(b)
                monoms += poly.monomials()
            monoms_list = tuple(Set(monoms))

            # check if the basis represented in terms of Monomials is efficient
            degs = [self._poly_ring(m).degree() for m in monoms]
            min_deg, max_deg = min(degs), max(degs)
            monoms_obj = Monomials(self._poly_ring, (min_deg, max_deg + 1))

            if monoms_obj.cardinality() < 2 * len(monoms_list):
                computational_basis = monoms_obj
            else:
                computational_basis = monoms_list
            self._monomial_module = PolynomialFreeModule(
                P=self._poly_ring, basis=computational_basis)
            cols = [self._monomial_module(b).to_vector() for b in self._basis]
            self._basis_mat = Matrix(cols).transpose()
            if self._basis_mat.ncols() > self._basis_mat.rank():
                raise ValueError(
                    "Basis polynomials are not linearly independent")

        def convert(self, p):
            r"""
            Algorithm is to convert all polynomials into monomials and use
            linear algebra to solve for the appropriate coefficients in this
            common basis.
            """
            try:
                p_vect = self._monomial_module(p).to_vector()
                decomp = self._basis_mat.solve_right(p_vect)
            except ValueError:
                raise ValueError(
                    "Value %s is not spanned by the basis polynomials" % p)
            polys = [v[1] * self._module.monomial(v[0])
                     for v in zip(self._basis, decomp)]
            module_p = sum(polys, self._module.zero())
            return module_p

    class _InfiniteBasisConverter:
        def __init__(self, P, comb_mod, basis):
            self._poly_ring = P
            self._module = comb_mod
            self._basis = basis
            # TODO implement some sort of caching for intermediate objects for
            # conversion

        def convert(self, p):
            # TODO check if polynomial is contained in span of basis in general
            # for now assume that infinite basis is in increasing order of
            # degree, and collect basis elements until the degree exceeds that
            # of input polynomial.
            # NOTE this is not necessarily a good heuristic without further
            # processing of the basis beforehand; always works for homogeneous
            # polynomials though
            # TODO for now, what characteristics of the basis guarantee
            # correctness here?
            deg = p.degree()
            fin_basis = []
            it = iter(self._basis)
            b = it.next()
            while b.degree() <= deg:
                fin_basis.append(b)
                b = it.next()

            C = PolynomialFreeModule._FiniteBasisConverter(
                self._poly_ring, self._module, fin_basis)
            return C.convert(p)

    def _element_constructor_(self, elt):
        return self.convert(elt)

    def convert(self, p):
        if p not in self._poly_ring:
            raise ValueError("Value %s is not a polynomial in %s"
                             % (p, self._poly_ring))
        # ensure that the type of p is actually a polynomial
        p = self._poly_ring(p)
        return self._converter.convert(p)

    Element = PolynomialFreeModuleElement

# TODO allow conversions between elements of different modules over the same
# polynomial space, and implement conversions/coercions

# would like to be able to:
# * convert between spaces seamlessly, find join and meet of spaces,
# * represent a basis of a polynomial space (as in an ideal with generators)
# * transition easily between representations in different bases

# perhaps it can automatically transition dimension, e.g. cast things up in
# dimension automatically but trim dimension only by explicit request
