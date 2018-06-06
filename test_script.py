from sage.all import *
import poly_utils
from zonotopal_algebra import ZonotopalAlgebra


def main():
    check_external_D_space(p401_cols, "xy")


simple_cols = [
    [1, 0],
    [0, 1],
    [1, 1]
]

A2_root_system_cols = [
    [1, -1,  0],
    [1,  0, -1],
    [0,  1, -1]
]

p401_cols = [
    [1, 0],
    [0, 1],
    [1, 1],
    [1, 0],
    [0, 1]
]

A3_root_system_cols = [
    [1, -1, 0, 0],
    [1, 0, -1, 0],
    [0, 1, -1, 0],
    [1, 0, 0, -1],
    [0, 1, 0, -1],
    [0, 0, 1, -1]
]

external_A3_root_system_cols = [
    [1, -1, 0, 0],
    [1, 0, -1, 0],
    [0, 1, -1, 0],
    [1, 0, 0, -1],
    [0, 1, 0, -1],
    [0, 0, 1, -1],
    [1, 0, 0, -1],
    [0, 1, 0, -1],
    [0, 0, 1, -1]
]

large_matroid_cols = [
    [1, -1, 0, 0],
    [1, 0, -1, 0],
    [0, 1, -1, 0],
    [1, 0, 0, -1],
    [0, 1, 0, -1],
    [0, 0, 1, -1],
    [1, 1, 1,  1],
    [1, 1, 1,  2],
    [1, 1, 1,  1],
    [2, 1, 3, -2],
    [2, 1, 3, -2]
]


def check_central_D_space(cols, varNames):
    X = Matrix(QQ, cols).transpose()
    Z = ZonotopalAlgebra(X, variant="central", data={'varNames': varNames})
    print "Z =", Z
    print
    print "Central P-space basis:"
    P_basis = Z.P_space_basis()
    for B in P_basis:
        print "%s: %s" % (tuple(B), P_basis[B].factor())
    print
    print "Central D-space basis:"
    D_basis = Z.D_space_basis()
    for B in D_basis:
        print "%s: %s" % (tuple(B), D_basis[B].factor())
    print
    print "Central I-ideal generators:"
    I_gens = Z.I_ideal_gens()
    for p in I_gens:
        print p.factor()
    print
    print "Central J-ideal generators:"
    J_gens = Z.J_ideal_gens()
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


def check_internal_D_space(cols, varNames):
    X = Matrix(QQ, cols).transpose()
    Z1 = ZonotopalAlgebra(X, variant="central", data={'varNames': varNames})
    Z2 = ZonotopalAlgebra(X, variant="internal", data={'varNames': varNames})
    print "Z2 =", Z2
    print
    print "Internal P-space basis:"
    P_basis = Z2.P_space_basis()
    for B in P_basis:
        print "%s: %s" % (tuple(B), P_basis[B].factor())
    print
    print "Internal D-space basis:"
    # central_D_basis = Z1.D_space_basis()
    # D_basis = {B: central_D_basis[B] for B in Z2._internal_bases()}
    D_basis = Z2.D_space_basis()
    for B in D_basis:
        print "%s: %s" % (tuple(B), D_basis[B].factor())
    print
    print "Internal I-ideal generators:"
    I_gens = Z2.I_ideal_gens()
    for p in I_gens:
        print p.factor()
    print
    print "Internal J-ideal generators:"
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

def check_external_D_space(cols, varNames):
    X = Matrix(QQ, cols).transpose()
    Z = ZonotopalAlgebra(X, variant="external", data={'varNames': varNames})
    print "Z =", Z
    print
    print "External P-space basis:"
    P_basis = Z.P_space_basis()
    for B in P_basis:
        print "%s: %s" % (tuple(B), P_basis[B].factor())
    print
    print "External D-space basis:"
    D_basis = Z.D_space_basis()
    for B in D_basis:
        print "%s: %s" % (tuple(B), D_basis[B].factor())
    print
    print "External I-ideal generators:"
    I_gens = Z.I_ideal_gens()
    for p in I_gens:
        print p.factor()
    print
    print "External J-ideal generators:"
    J_gens = Z.J_ideal_gens()
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


# def main():
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
