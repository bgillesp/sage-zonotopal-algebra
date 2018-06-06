from central_zonotopal_algebra import CentralZonotopalAlgebra
from external_zonotopal_algebra import ExternalZonotopalAlgebra
from internal_zonotopal_algebra import InternalZonotopalAlgebra
from sage.repl.rich_output.pretty_print import pretty_print


def ZonotopalAlgebra(X, variant="central", data={}):
    if variant == "central":
        return CentralZonotopalAlgebra(X, **data)
    elif variant == "external":
        return ExternalZonotopalAlgebra(X, **data)
    elif variant == "internal":
        return InternalZonotopalAlgebra(X, **data)
    else:
        raise ValueError("uncrecognized zonotopal algebra type: %s" % variant)


def zon_spaces(Z, spaces="IJPD"):
    I, J, P, D = None, None, None, None
    if "I" in spaces:
        print "Generating I ideal gens..."
        I = Z.I_ideal_gens()
        I.sort()
        I = [i.factor() for i in I]
    if "J" in spaces:
        print "Generating J ideal gens..."
        J = Z.J_ideal_gens()
        J.sort()
        J = [j.factor() for j in J]
    if "P" in spaces:
        print "Generating P space basis..."
        P = Z.P_basis()
        P.sort()
        P = [p.factor() for p in P]
    if "D" in spaces:
        print "Generating D space basis..."
        D = Z.D_basis()
        D.sort()
        D = [d.factor() for d in D]
    return (I, J, P, D)


def print_zon_info(tup):
    I, J, P, D = tup
    if I is not None:
        print "I(X) ="
        pretty_print(I)
        print
    if J is not None:
        print "J(X) ="
        pretty_print(J)
        print
    if P is not None:
        print "P(X) ="
        pretty_print(P)
        print
    if D is not None:
        print "D(X) ="
        pretty_print(D)
        print


def zon_info(X, variant="central", data={}, spaces="IJPD"):
    data = data.copy()
    if 'varNames' not in data:
        rws = X.nrows()
        varstr = "xyzwabcd"
        if rws <= len(varstr):
            varstr = varstr[:rws]
        else:
            varstr = varstr[0]
        data['varNames'] = varstr
    Z = ZonotopalAlgebra(X, variant, data)
    tup = zon_spaces(Z, spaces)
    # print out central zonotopal spaces
    print "X = "
    pretty_print(Z.matrix())
    print ""
    print_zon_info(tup)
    return tup
