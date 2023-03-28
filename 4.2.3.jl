using SymPy

yVal = .93
@vars X
expr = X^2 / (1 + X^2) - yVal
soln = solve(expr,X)
##COMPARE TO THE OTHER RESULTS