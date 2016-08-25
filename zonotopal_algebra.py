class ZonotopalAlgebra:
    def __init__(self,X,variant="central",data={}):
        self.variant = variant
        self.data = data
        self.obj = None
        if type=="central":
            varNames = 'x'; # TODO default variable
            if 'varNames' in data:
                varNames = data['varNames']
            self.obj = CentralZonotopalAlgebra(X, varNames)
        # TODO include additional types of zonotopal algebras
        else:
            raise ValueException("uncrecognized zonotopal algebra type: %s" % type)

    def I_ideal_gens(self):
        return self.obj.I_ideal_gens()

    def J_ideal_gens(self):
        return self.obj.I_ideal_gens()

    def D_basis(self):
        return self.obj.D_basis()

    def P_basis(self):
        return self.obj.P_basis()

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
