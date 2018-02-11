from central_zonotopal_algebra import CentralZonotopalAlgebra
from internal_zonotopal_algebra import InternalZonotopalAlgebra
from external_zonotopal_algebra import ExternalZonotopalAlgebra
from forward_exchange_zonotopal_algebra import ForwardExchangeZonotopalAlgebra
from sage.repl.rich_output.pretty_print import pretty_print

# TODO could reformulate this as a factory method rather than a shim
def ZonotopalAlgebra(X,variant="central",data={}):
    if variant == "central":
        return CentralZonotopalAlgebra(X, **data)
    elif variant == "internal":
        return InternalZonotopalAlgebra(X, **data)
    elif variant == "external":
        return ExternalZonotopalAlgebra(X, **data)
    elif variant == "forward exchange":
        return ForwardExchangeZonotopalAlgebra(X, **data)
    # TODO include additional types of zonotopal algebras as possible
    else:
        raise ValueError("uncrecognized zonotopal algebra type: %s" % variant)

    # def __repr__(self):
    #     return self.obj.__repr__()
    #
    # def I_ideal_gens(self):
    #     return self.obj.I_ideal_gens()
    #
    # def I_ideal(self):
    #     return self.obj.I_ideal()
    #
    # def J_ideal_gens(self):
    #     return self.obj.J_ideal_gens()
    #
    # def J_ideal(self):
    #     return self.obj.J_ideal()
    #
    # def D_basis(self):
    #     return self.obj.D_basis()
    #
    # def P_basis(self):
    #     return self.obj.P_basis()

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

def zon_info(X, variant="central", data={}):
    data = data.copy()
    if not varNames in data:
        rws = X.nrows()
        varstr = "xyzabc"
        if rws <= len(varstr):
            varstr = varstr[:rws]
        else:
            varstr = varstr[0]
        data.varNames = varstr
    Z = ZonotopalAlgebra(X, variant, data)
    tup = zon_spaces(Z)
    # print out central zonotopal spaces
    print "X = "
    pretty_print(Z.matrix())
    print ""
    print_zon_info(tup)
    return tup
