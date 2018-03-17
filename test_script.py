from sage.all import *

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
    #
    # cols2 = [
    #     [1, 0],
    #     [1, 1],
    #     [0, 1],
    #     [0, 1]
    # ]
    # X2 = Matrix(QQ, cols2).transpose()
    #
    # Z2 = ZonotopalAlgebra(X2, variant="internal")
    # print "Z2 =", Z2
    # print "Internal P space basis:"
    # print Z2.P_space_basis()

    cols3 = [
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 1],
        [1, 0, 0],
        [1, 1, 0]
    ]
    X3 = Matrix(QQ, cols3).transpose()

    Z3 = ZonotopalAlgebra(X3, variant="internal")
    print "Z3 =", Z3
    print "Internal bases:"
    print list(Z3._internal_bases())
    print "Hyperplanes"
    print list(Z3._ordered_matroid().hyperplanes())
    gens = Z3.I_ideal_gens()
    for g in gens:
        print factor(g)


if __name__ == "__main__":
    main()
