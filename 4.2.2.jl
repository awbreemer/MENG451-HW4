using SymPy

x = [0:5;]
y = [0 .5 .8 .9 .941176 .961538]'

xForThis = x[2:5]
yForthis = y[2:5]

function lagInterp(x, y)

    @vars X

    ls = ones(typeof(X), length(x))
    for j in eachindex(x)
        for m in eachindex(x)
            if j != m
                ls[j] *= (X - x[m]) / (x[j] - x[m])
            end
        end
    end

    L :: typeof(X)= 0 
    for k in eachindex(x)
        L += y[k] * ls[k]
    end

    L = simplify(L)
    pol = sympy.Poly(L, X)
    coefs = pol.coeffs()
    return coefs, L
end

lagPoly = lagInterp(xForThis, yForthis)

length(lagPoly[1])

cofs  = zeros(Float64, 4)

for i in eachindex(lagPoly[1])
    cofs[i] = lagPoly[1][i] 
end

function solveEqn(eqn, equals)
    @vars X
    eqn -= equals
    ans = solve(eqn, X)
    return ans
end


#Using the SymPy Equation solver
xAnswer = solveEqn(lagPoly[2], .93)

##Using newtons method
