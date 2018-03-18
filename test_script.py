from sage.all import *

import poly_utils
from zonotopal_algebra import ZonotopalAlgebra


def main():
    # cols1 = [
    #     [1, -1, 0, 0],
    #     [1, 0, -1, 0],
    #     [0, 1, -1, 0],
    #     [1, 0, 0, -1],
    #     [0, 1, 0, -1],
    #     [0, 0, 1, -1]
    # ]
    # X1 = Matrix(QQ, cols1).transpose()
    #
    # Z1 = ZonotopalAlgebra(X1, variant="internal")
    # P1 = Z1.polynomial_ring()
    # x0, x1, x2, x3 = P1.gens()
    # print "Z1 =", Z1
    # p = (x1 - 2*x2 + x3)
    # print "Polynomial:", p
    # print "In P-space basis:", Z1.P_space()(p)

    # cols2 = [
    #     [1, 0],
    #     [0, 1],
    #     [1, 1],
    #     [1, 0],
    #     [0, 1]
    # ]
    cols2 = [
        [1, -1, 0, 0],
        [1, 0, -1, 0],
        [0, 1, -1, 0],
        [1, 0, 0, -1],
        [0, 1, 0, -1],
        [0, 0, 1, -1],
        [1, 1, 1,  1],
        [1, 1, 1,  2],
        [1, 1, 1,  1]
    ]
    X2 = Matrix(QQ, cols2).transpose()

    Z2 = ZonotopalAlgebra(X2, variant="central", data={'varNames': "xyzw"})
    print "Z2 =", Z2
    print
    print "Central P-space basis:"
    P_basis = Z2.P_space_basis()
    for B in P_basis:
        print "%s: %s" % (tuple(B), P_basis[B].factor())
    print
    print "Central D-space basis:"
    D_basis = Z2.D_space_basis()
    for B in D_basis:
        print "%s: %s" % (tuple(B), D_basis[B].factor())
    print
    print "Central I-ideal generators:"
    I_gens = Z2.I_ideal_gens()
    for p in I_gens:
        print p.factor()
    print
    print "Central J-ideal generators:"
    J_gens = Z2.J_ideal_gens()
    for p in J_gens:
        print p.factor()
    print
    print "P-D duality check:"
    for B1 in P_basis:
        p = P_basis[B1]
        for B2 in D_basis:
            d = D_basis[B2]
            print poly_utils.diff_bilinear_form(p, d),
        print
    print
    print "D = ker(J) check:"
    for d in D_basis.values():
        for j in J_gens:
            print poly_utils.poly_deriv(j, d),
        print

    # cols3 = [
    #     [0, 0, 1],
    #     [0, 1, 0],
    #     [1, 0, 1],
    #     [1, 0, 0],
    #     [1, 1, 0]
    # ]
    # X3 = Matrix(QQ, cols3).transpose()
    #
    # Z3 = ZonotopalAlgebra(X3, variant="internal")
    # print "Z3 =", Z3
    # print "Internal bases:"
    # print list(Z3._internal_bases())
    # print "Hyperplanes"
    # print list(Z3._ordered_matroid().hyperplanes())
    # gens = Z3.I_ideal_gens()
    # for g in gens:
    #     print factor(g)


if __name__ == "__main__":
    main()
