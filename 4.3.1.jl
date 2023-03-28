using SymPy
using Polynomials

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

f(x) = 1 / (1 + 25 * x^2)

N = 11 #Number of evenly spaced points
x = LinRange(-1, 1, N) #x value points
y = f.(x) #y value points

lagPoly = lagInterp(x, y)

cofs = zeros(Float64, length(lagPoly[1]))
for i in eachindex(lagPoly[1])
    cofs[i] = lagPoly[1][i] 
end

thePoly = Polynomial(reverse(cofs))

resultXs = LinRange(-1,1, 100)
thePolyVals = thePoly.(resultXs)
theActualVals = f.(resultXs)

t = LinRange(0, pi, N)
chebyshevNodes = -cos.(t)
chebyYs = f.(chebyshevNodes)
chebyLag = lagInterp(chebyshevNodes, chebyYs)

cofsC = zeros(Float64, length(chebyLag[1]))
for i in eachindex(chebyLag[1])
    cofsC[i] = chebyLag[1][i] 
end

chebyPoly = Polynomial(reverse(cofsC))
chebyVals = chebyPoly.(resultXs)

plot(resultXs, [theActualVals, thePolyVals, chebyVals])