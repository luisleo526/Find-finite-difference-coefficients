import numpy as np
from math import floor, factorial
from sympy.matrices import Matrix
from sympy import Symbol, Add, Eq
from sympy import init_printing, pprint, latex
from sympy.printing.str import StrPrinter

def f(i):
    if i == 0:
        return Symbol('f_i')
    else:
        return Symbol('f_{i%+d}' % i)


def central_coeff(m, n):
    assert n % 2 == 0, "Wrong order"
    p = int((2 * floor((m + 1) / 2) - 1 + n - 1) / 2)
    A = np.zeros((2 * p + 1, 2 * p + 1), dtype=int)
    for i in range(2 * p + 1):
        for j in range(2 * p + 1):
            A[i, j] = (-p + j) ** i
    b = np.zeros(2 * p + 1, dtype=int)
    b[m] = factorial(m)
    A = Matrix(A)
    b = Matrix(b)
    x = A.LUsolve(b)
    expr = []
    for i in range(2 * p + 1):
        expr.append(f(-p + i) * x[i])
    expr = StrPrinter(dict(order='none'))._print_Add(Add(*expr, evaluate=False))
    eq = Eq(Symbol('d^{%d} f' % m) / Symbol('d x^{%d}' % m), Symbol(expr) / Symbol('h^%d' % m),
            evaluate=False)
    return eq


def side(s, d):
    n = len(s)
    assert n > d
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            A[i, j] = s[j] ** i
    b = np.zeros(n, dtype=int)
    b[d] = factorial(d)
    A = Matrix(A)
    b = Matrix(b)
    x = A.LUsolve(b)
    expr = []
    for i in range(n):
        expr.append(f(s[i]) * x[i])
    expr = StrPrinter(dict(order='none'))._print_Add(Add(*expr, evaluate=False))
    eq = Eq(Symbol('d^%d f' % d) / Symbol('d x^%d' % d), Symbol(expr) / Symbol('h^%d' % d),
            evaluate=False)
    return eq

if __name__ == '__main__':

    init_printing()

    method = int(input("Choose forward/backward (0) or central (1) method: "))
    assert method in [0, 1], "Invalid input"
    order = int(input("Input order of derivatives: "))
    assert order > 0, "Order should be positive"
    if method == 0:
        indices = [int(x) for x in input("Input indices, separate with ','").strip().split(',')]
        print("")
        pprint(side(indices, order), use_unicode=False)
    else:
        n=int(input("Order accuracy: "))
        print("")
        pprint(central_coeff(order, n), use_unicode=False)
