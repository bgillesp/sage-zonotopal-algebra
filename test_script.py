from sage.all import *

import zonotopal_algebra
from ordered_matroid import OrderedMatroid


def main():
    cols = [[1, 0], [0, 1], [1, 1], [1, 1]]
    X = Matrix(QQ, cols).transpose()
    # M = OrderedMatroid(Matroid(X), reverse=True)
    # print M
    # print M.coactive_elements([2])
    zonotopal_algebra.zon_info(X, "internal", spaces="IJP")


if __name__ == "__main__":
    main()
