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
        return "Zonotopal Algebra over the field " + str(self.F) + " with matrix " + str(self.X)

    def _poly(self,v):
        n = self.Pi.ngens()
        vars = self.Pi.gens()
        p = self.Pi.zero()
        for i in range(n):
            p += v[i]*vars[i]
        return p

    def _polys(self,S,matr):
        p = self.Pi.one()
        for s in S:
            p *= self._poly(matr.column(s))
        return p

    def __facet_hyperplanes(self):
        if self.hyperplanes == None:
            subsets = []
            eta = []
            hyp = self.M.hyperplanes() # operation gets slow for large X
            for h in hyp:
                subsets.append(h)
                vectors = self.M.max_independent(h)
                vectors = map(lambda v: X.column(v), vectors)
                space = (self.V.subspace(vectors)).complement()
                eta.append(space.basis()[0])
            self.hyperplanes = (eta, subsets) # zippable
        return self.hyperplanes

    def _facet_hyperplanes(self):
        return self.__facet_hyperplanes()[1]

    def _facet_hyperplane_normals(self):
        return self.__facet_hyperplanes()[0]

    def I_ideal_gens(self):
        if self.I == None:
            gens = []
            eta = self._facet_hyperplane_normals()
            subsets = self._facet_hyperplanes()
            for i in range(len(eta)):
                gens.append(self._poly(eta[i])**(self.M.size()-len(subsets[i])))
            self.I = gens
        return self.I

    def I_ideal(self):
        return ideal(self.I_ideal_gens())

    # not used currently for computation of zonotopal spaces
    def _short_long_sets(self):
        if self.short_long == None:
            S = []; L = []
            # probably can do this more efficiently
            for s in powerset(self.M.groundset()):
                if self.M.is_coindependent(s):
                    S.append(s)
                else:
                    L.append(s)
            self.short_long = (S, L)
        return self.short_long

    def short_subsets(self):
        return Sequence(self._short_long_sets()[0], immutable=True)

    def long_subsets(self):
        return Sequence(self._short_long_sets()[1], immutable=True)

    def _cocircuits(self):
        if self.cocircuits == None:
            # complements of hyperplanes, this construction is zippable with corresponding hyperplanes
            cocircs = []
            hyperplanes = self._facet_hyperplanes()
            g = self.M.groundset()
            for h in hyperplanes:
                cocircs.append(g.difference(h))
            self.cocircuits = cocircs
        return self.cocircuits

    def J_ideal_gens(self):
        if self.J == None:
            gens = []
            for c in self._cocircuits():
                gens.append(self._polys(c,self.X))
            self.J = gens
        return self.J

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

    def P_basis(self):
        if self.P == None:
            P = []
            for b in self.M.bases():
                elts = self._big_ex(b)
                P.append(self._polys(elts,self.X))
            P.sort()
            self.P = P
        return self.P

    def _poly_dual_basis(self, poly_basis):
        # compute a dual basis for an input *homogeneous* basis
        class _PolynomialVectorSpace(CombinatorialFreeModule):
            def __init__(self, R, n, k):
                self._name = 'Module of degree at most %d polynomials in %d variables' % (n,k)
                CombinatorialFreeModule.__init__(self, R, IntegerListsLex(max_sum=n, length=k))
                self.print_options(prefix='x', bracket=['^', ''], latex_bracket=['^{', '}'])
        deg = max([p.degree() for p in poly_basis])
        poly_vectors = _PolynomialVectorSpace(self.F, deg, self.X.nrows())
        converted_basis = []
        for p in poly_basis:
            poly_items = [(poly_vectors(exponent), coeff) for exponent, coeff in p.dict().items()]
            v = poly_vectors.linear_combination(poly_items)
            converted_basis.append(v)
        # compute dual basis
        bilinear_form_coeffs = \
            [prod(map(factorial, b.monomial_coefficients().keys()[0])) for b in poly_vectors.basis()]
        A = Matrix([v.to_vector() for v in converted_basis])
        D = Matrix.diagonal(bilinear_form_coeffs)
        B = (A*D*A.transpose()).inverse()
        dual_poly_basis = []
        for col in B.columns():
            q = sum(coeff*p for coeff,p in zip(col, poly_basis))
            dual_poly_basis.append(q)
        return dual_poly_basis

    def _D_basis_project(self,dual_basis):
        D = []
        # associate hyperplanes with cocircuits
        hyperplane_normals = self._facet_hyperplane_normals()
        cocircuits = self._cocircuits()
        polys = [] # hyperplane normals and cocircuits are ordered to allow zipping
        for eta, c in zip(hyperplane_normals, cocircuits):
            i = self._poly(eta)**len(c)
            j = self._polys(c,self.X)
            polys.append( (i,j) )
        for p in dual_basis:
            q = p
            for (i,j) in polys:
                coeff1 = CentralZonotopalAlgebra.diff_bilinear_form(j,p)
                coeff2 = CentralZonotopalAlgebra.diff_bilinear_form(j,i)
                q -= coeff1/coeff2 * i
            D.append(q)
        return D

    def D_basis(self):
        # modify P(X) basis (acting as a P(X)* basis) to be in kernel of J(X)
        if self.D == None:
            P_basis = self.P_basis()
            P_basis_internal_dual = self._poly_dual_basis(P_basis)
            self.D = self._D_basis_project(P_basis_internal_dual)
        return self.D

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
        n_vars = len(p.parent().gens())
        zero = [0]*n_vars
        return (poly_deriv(p,q))(zero)

def zon_spaces(Z):
    print "Generating I..."
    #I = Sequence(Z.I_ideal().gens(),universe=Z.Pi)
    I = Z.I_ideal_gens()
    I.sort()
    print "Generating J..."
    #J = Sequence(Z.J_ideal().groebner_basis(),universe=Z.Pi)
    J = Z.J_ideal_gens()
    J.sort()
    print "Generating P..."
    P = Z.P_basis()
    P.sort()
    print "Generating D..."
    D = Z.D_basis()
    D.sort()
    return (I,J,P,D)

def print_zon_info(tup):
    (I,J,P,D) = tup
    print "I(X) ="
    pretty_print(I)
    print ""
    print "J(X) ="
    pretty_print(J)
    print ""
    print "P(X) ="
    pretty_print(P)
    print ""
    print "D(X) ="
    pretty_print(D)

def zon_info(X):
    rws = X.nrows()
    varstr = "xyzabc"
    if rws <= len(varstr):
        varstr = varstr[:rws]
    else:
        varstr = varstr[0]
    Z = CentralZonotopalAlgebra(X,varstr)
    tup = zon_spaces(Z)
    # print out central zonotopal spaces
    print "X = "
    pretty_print(Z.X)
    print ""
    print_zon_info(tup)
    return tup
