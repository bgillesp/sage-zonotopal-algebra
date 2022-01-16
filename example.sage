from zonotopal_algebra import *

U23_cols = [
    [1, 0],
    [0, 1],
    [1, 1],
]
U23 = Matrix(QQ, U23_cols).transpose()

Z1 = ZonotopalAlgebra(
    U23,
    variant="central",
    varNames="xy"
)
print("Central zonotopal spaces of U^2_3:")
print()
print_zon_info(zon_spaces(Z1))

Z2 = ZonotopalAlgebra(
    U23,
    variant="external",
    varNames="xy"
)
print("External zonotopal spaces of U^2_3:")
print()
print_zon_info(zon_spaces(Z2))

Z3 = ZonotopalAlgebra(
    U23,
    variant="internal",
    varNames="xy"
)
print("Internal zonotopal spaces of U^2_3:")
print()
print_zon_info(zon_spaces(Z3))
