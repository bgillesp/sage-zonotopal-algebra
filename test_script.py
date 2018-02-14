from sage.all import *

import zonotopal_algebra
from ordered_matroid import OrderedMatroid


def main():
    cols = [[1, 0], [0, 1], [1, 1]]
    X = Matrix(QQ, cols).transpose()
    # M = OrderedMatroid(Matroid(X), reverse=True)
    # print M
    # print M.coactive_elements([2])
    zonotopal_algebra.zon_info(X, "external", spaces="IJPD")
    # Z = zonotopal_algebra.ZonotopalAlgebra(X, "external")
    # print "Z = ", Z
    # print "External Matrix =\n", Z._extending_basis_matrix()
    # print Z._external_cocircuit([2,3])
    # print Z._external_basis([1,2])


if __name__ == "__main__":
    main()
