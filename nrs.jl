"""nrsδ(x,f,f′)
calculate the next correction to x, our estimate of a root of f, with
derivative f′
"""
nrsδ(x, f, f′) = - f(x) / f′(x)

function nrs(x_0::Number, f, f′)
    ϵ = 10^(-6)
    count = 0
    x_n = x_0
    Δ = 1.0
    while Δ > ϵ && count < 40
        Δ = nrsδ(x_n, f, f′)
        x_n = x_n + Δ
        count = count + 1
    end
    return x_n
end


F(x) = x^3 - 8
dF(x) = 3x^2

println(nrs(5, F, dF))


function nrs(xs::Vector{T}, f, f′) where T <: Number
    """
    Take a numeric vector as first argument instead,
    to preserve history of each iteration.
    """
    ϵ = 10^(-6)
    count = 0
    x_n = xs[1]
    Δ = 1.0
    while abs(Δ) > ϵ && count < 40
        Δ = nrsδ(x_n, f, f′)
        x_n = x_n + Δ
        push!(xs, x_n)
        count = count + 1
    end
    return xs
end

println(nrs([5.0], F, dF))

function nrs(xs::Vector{T}, f, f′) where T <: Number
    """
    nrs which does not modify its Vector argument
    """
    xsc = deepcopy(xs)
    ϵ = 10^(-6)
    count = 0
    x_n = xs[1]
    Δ = 1.0
    while abs(Δ) > ϵ && count < 40
        Δ = nrsδ(x_n, f, f′)
        x_n = x_n + Δ
        push!(xsc, x_n)
        count = count + 1
    end
    return xsc
end

xsa = [5.0]
println(nrs.(xsa, F, dF))
println(xsa)

using Plots

Plots.plot(0.0:0.01:3.0, F)

history = nrs([1.2], F, dF)

println(history)
println(F.(history))

# plot! adds a plot to the existing plotting window
Plots.plot!(history, F, seriestype=:scatter)

z(x, y) = abs(F(x + im*y))
Plots.contour(-3.0:0.01:3.0, -3.0:0.01:3.0, z)

numb = complex([5.0 + im*2])
println(nrs(numb, F, dF))
historyC = nrs([1.0+im*3], F, dF)
println(real.(historyC))
println(imag.(historyC))
Plots.plot!(real.(historyC), imag.(historyC), seriestype=:scatter, markercolor=1:10)
