from polynomial_vector_space import PolynomialModule
from polynomial_vector_space import Monomials
from poly_utils import PolyUtils
from sage.misc.cachefunc import cached_method
from sage.matroids.constructor import Matroid
from sage.modules.free_module import VectorSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod
from sage.functions.other import factorial
from sage.matrix.constructor import Matrix

class CentralZonotopalAlgebra(object):
    def __init__(self,X,varNames="x"):
        self.F = X.base_ring()
        self.M = Matroid(matrix=X)
        self.X = self.M.representation()
        self.V = VectorSpace(self.F,self.X.nrows())
        self.Pi = PolynomialRing(self.F, self.X.nrows(), names=varNames, order='deglex')
        # these are computed as-needed
        self.hyperplanes = None # hyperplanes, as a tuple of lists, (orthogonal vectors, hyperplane subsets)
        self.short_long = None # short and long subsets
        self.cocircuits = None
        self.I = None # I(X) ideal generators
        self.J = None # J(X) ideal generators
        self.P = None # P(X) homogeneous basis elements
        self.D = None # D(X) homogeneous basis elements

    def __repr__(self):
        return "Zonotopal Algebra over the field " + str(self.F) + " with matrix\n" + str(self.X)

    @cached_method
    def __facet_hyperplanes(self):
        subsets = []
        eta = []
        hyp = self.M.hyperplanes() # operation gets slow for large X
        for h in hyp:
            subsets.append(h)
            vectors = self.M.max_independent(h)
            vectors = map(lambda v: self.X.column(v), vectors)
            space = (self.V.subspace(vectors)).complement()
            eta.append(space.basis()[0])
        return (eta, subsets) # zippable

    def _facet_hyperplanes(self):
        return self.__facet_hyperplanes()[1]

    def _facet_hyperplane_normals(self):
        return self.__facet_hyperplanes()[0]

    @cached_method
    def I_ideal_gens(self):
        gens = []
        eta = self._facet_hyperplane_normals()
        subsets = self._facet_hyperplanes()
        for i in range(len(eta)):
            gen = PolyUtils.linear_form(self.Pi, eta[i])**(self.M.size() - len(subsets[i]))
            gens.append(gen)
        return gens

    def I_ideal(self):
        return ideal(self.I_ideal_gens())

    # not used currently for computation of zonotopal spaces
    @cached_method
    def _short_long_sets(self):
        S = []; L = []
        # probably can do this more efficiently
        for s in powerset(self.M.groundset()):
            if self.M.is_coindependent(s):
                S.append(s)
            else:
                L.append(s)
        return (S, L)

    def short_subsets(self):
        return Sequence(self._short_long_sets()[0], immutable=True)

    def long_subsets(self):
        return Sequence(self._short_long_sets()[1], immutable=True)

    @cached_method
    def _cocircuits(self):
        # complements of hyperplanes, this construction is zippable with corresponding hyperplanes
        cocircs = []
        hyperplanes = self._facet_hyperplanes()
        g = self.M.groundset()
        for h in hyperplanes:
            cocircs.append(g.difference(h))
        return cocircs

    @cached_method
    def J_ideal_gens(self):
        gens = []
        for c in self._cocircuits():
            gen = PolyUtils.pure_tensor(self.Pi, self.X, c)
            gens.append(gen)
        return gens

    def J_ideal(self):
        return ideal(self.J_ideal_gens())

    def _big_ex(self,S):
        ex_poss = self.M.groundset().difference(S)
        def f(e):
            flat = self.M.closure(filter(lambda i: i <= e, S))
            return not e in flat
        ex = filter(f, ex_poss)
        return ex

    def _big_y(self,S):
        y_poss = self.M.groundset().difference(S)
        def f(e):
            flat = self.M.closure(filter(lambda i: i >= e, S))
            return not e in flat
        y = filter(f, y_poss)
        return y

    @cached_method
    def P_basis(self):
        P = []
        for b in self.M.bases():
            xb = self._big_ex(b)
            elt = PolyUtils.pure_tensor(self.Pi, self.X, xb)
            P.append(elt)
        P.sort()
        return P

    def _poly_dual_basis(self, poly_basis):
        # compute a dual basis for an input *homogeneous* basis
        deg = max([p.degree() for p in poly_basis])
        poly_module = PolynomialModule( self.Pi, basis=Monomials(self.Pi, (0, deg+1)) )
        bilinear_form_coeffs = []
        for b in poly_module.basis().keys():
            # each b is a monomial in self.Pi of degree at most deg
            b = self.Pi(b)
            bilinear_form_coeffs.append( prod( map(factorial, b.degrees()) ) )
        A = Matrix([poly_module(p).to_vector() for p in poly_basis])
        D = Matrix.diagonal(bilinear_form_coeffs, sparse=False)
        B = (A*D*A.transpose()).inverse()
        dual_basis = []
        for col in B.columns():
            q = sum( [coeff*p for coeff,p in zip(col, poly_basis)] )
            dual_basis.append(q)
        return dual_basis

    def _D_basis_project(self,dual_basis):
        D = []
        # associate hyperplanes with cocircuits
        hyperplane_normals = self._facet_hyperplane_normals()
        cocircuits = self._cocircuits()
        polys = [] # hyperplane normals and cocircuits are ordered to allow zipping
        for eta, c in zip(hyperplane_normals, cocircuits):
            i = PolyUtils.linear_form(self.Pi, eta)**len(c)
            j = PolyUtils.pure_tensor(self.Pi, self.X, c)
            polys.append( (i,j) )
        for p in dual_basis:
            q = p
            for (i,j) in polys:
                coeff1 = PolyUtils.diff_bilinear_form(j,p)
                coeff2 = PolyUtils.diff_bilinear_form(j,i)
                q -= coeff1/coeff2 * i
            D.append(q)
        return D

    @cached_method
    def D_basis(self):
        # modify P(X) basis (acting as a P(X)* basis) to be in kernel of J(X)
        P_basis = self.P_basis()
        P_basis_internal_dual = self._poly_dual_basis(P_basis)
        return self._D_basis_project(P_basis_internal_dual)
