from monomials import Monomials
from poly_free_module import PolynomialFreeModule
from sage.functions.other import factorial
from sage.matrix.constructor import Matrix
from sage.misc.misc_c import prod


def poly_deriv(p, q):
    r"""
    The derivative of q by the polynomial p, treating variables in p as
    corresponding partial differential operators.

    INPUT:

    - ``p`` -- a polynomial to act as a differential operator
    - ``q`` -- a polynomial to be differentiated

    OUTPUT:

    - the polynomial derivative of ``q`` by the operator corresponding to ``p``

    EXAMPLES:

        sage: P.<x, y> = PolynomialRing(QQ)
        sage: poly_deriv(x, x)
        1
        sage: poly_deriv(P.one(), x^4 + 2*x*y + 1)
        x^4 + 2*x*y + 1
        sage: poly_deriv(y, x)
        0
        sage: poly_deriv(x+y, x-y)
        0
        sage: poly_deriv(3*x^2, x^3 + x^2*y + x^2)
        18*x + 6*y + 6
    """
    g = p.parent().gens()
    s = p.parent().zero()
    for e_tup, coeff in p.dict().iteritems():
        diff_list = []
        for v, e in zip(g, e_tup):
            diff_list.extend([v] * e)
        s += coeff * q.derivative(diff_list)
    return s


def diff_bilinear_form(p, q):
    """
    Return the differential bilinear form `<p|q>` of ``p`` with ``q``

    The polynomial differential bilinear form is given by differentiating
    ``q`` by ``p`` as in the function ``poly_deriv``, and then evaluating the
    result at `x = 0`.

    Although not immediate from the definition, the operation is symmetric with
    respect to ``p`` and ``q``.

    The monomial basis is orthogonal with respect to this bilinear form, and
    the squared norm of a monomial is given by the product of factorials of its
    exponents.

    INPUT:

    - ``p`` -- a polynomial
    - ``q`` -- a polynomial

    OUTPUT:

    - the differential bilinear form of p and q

    EXAMPLES:

        sage: P.<x, y> = PolynomialRing(QQ)
        sage: diff_bilinear_form(x^2, x^2)
        2
        sage: diff_bilinear_form(x^3*y, x^3*y)
        6
        sage: diff_bilinear_form(x + y^2, x + 2*y + y^2)
        3
    """
    deriv_poly = poly_deriv(p, q)
    n_vars = len(p.parent().gens())
    zero = [0] * n_vars
    return deriv_poly(zero)


def linear_form(P, vec):
    r"""
    Return a linear form in the polynomial ring

    INPUT:

    - ``P`` -- a polynomial ring
    - ``vec`` -- a list of coefficients

    OUTPUT:

    - the linear form in ``P`` whose coefficients are given by ``vec``, where
      the ordering on the coefficients corresponds with the ordering of
      polynomial generators given by ``P.gens()``

    EXAMPLE:

        sage: P.<x, y> = PolynomialRing(QQ)
        sage: linear_form(P, [1, 2])
        x + 2*y
    """
    terms = [coeff * var for coeff, var in zip(vec, P.gens())]
    return sum(terms, P.zero())


def pure_tensor(P, vects, indices=None):
    r"""
    Return the product of a collection of linear forms corresponding to vectors

    INPUT:

    - ``P`` -- a polynomial ring
    - ``vects`` -- a list of coefficient vectors
    - ``indices`` -- (default: ``None``) a collection of indices for the list
      ``vects`` corresponding to the linear forms which will be included in the
      product.  Indices may be repeated.  If ``None``, then each vector will be
      included in the product once.

    OUTPUT:

    - the corresponding product of linear forms in the polynomial ring ``P``

    EXAMPLES:

        sage: P.<x, y> = PolynomialRing(QQ)
        sage: vects = [[1,0],[0,1],[1,1]]
        sage: pure_tensor(P, vects).factor()
        y * x * (x + y)
        sage: pure_tensor(P, vects, [0, 2, 2]).factor()
        x * (x + y)^2
    """
    if indices is None:
        indices = range(len(vects))
    terms = [linear_form(P, vects[i]) for i in indices]
    return prod(terms, P.one())


def poly_dual_basis(P, poly_basis):
    r"""
    Return a collection of polynomials which are dual under the differential
    bilinear form to a given homogeneous collection

    INPUT:

    - ``P`` -- a polynomial ring
    - ``poly_basis`` -- a collection of polynomials in ``P`` which are
      homogeneous and linearly independent

    OUTPUT:

    - the dual basis of the polynomials in ``poly_basis`` in their span

    EXAMPLES:

        sage: P.<x, y> = PolynomialRing(QQ)
        sage: poly_basis = (1, x, x+y)
        sage: poly_dual_basis(P, poly_basis)
        [1, x - y, y]
        sage: poly_basis = (1, 2*x - y, x^2, x^2 + x*y)
        sage: poly_dual_basis(P, poly_basis)
        [1, 2/5*x - 1/5*y, 1/2*x^2 - x*y, x*y]
    """
    # recast poly_basis to ensure elements are all from P
    poly_basis = [P(p) for p in poly_basis]
    # compute max degree of basis polynomials for linear algebra computations
    deg = max([p.degree() for p in poly_basis])
    # construct polynomial free module for linear algebra computations
    monoms = Monomials(P, (0, deg))
    poly_module = PolynomialFreeModule(P, basis=monoms)
    # compute the values of the bilinear form <m|m> for basis monomials m
    bilinear_form_coeffs = []
    for b in poly_module.basis().keys():
        # each b is a monomial in P of degree at most deg
        b = P(b)
        bilinear_form_coeffs.append(prod(map(factorial, b.degrees())))
    # compute dual basis
    A = Matrix([poly_module(p).to_vector() for p in poly_basis])
    D = Matrix.diagonal(bilinear_form_coeffs, sparse=False)
    B = (A * D * A.transpose()).inverse()
    # reconstruct dual basis polynomials from corresponding vectors
    dual_basis = []
    for col in B.columns():
        q = sum([coeff * p for coeff, p in zip(col, poly_basis)])
        dual_basis.append(q)
    return dual_basis
