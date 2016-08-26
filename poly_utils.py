from sage.misc.misc_c import prod

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
        n_vars = len(p.parent().gens())
        zero = [0]*n_vars
        return (PolyUtils.poly_deriv(p,q))(zero)

    @staticmethod
    def linear_form(P, vec):
        terms = [coeff*var for coeff, var in zip(vec, P.gens())]
        return sum(terms, P.zero())

    @staticmethod
    def pure_tensor(P, matr, indices):
        terms = [PolyUtils.linear_form(P, matr.column(i)) for i in indices]
        return prod(terms, P.one())
