using Polynomials

include("plot.jl")

T = 0.5
K = 1.0

function polynomial_roots(coefficients)
    return roots(Polynomial(coefficients))
end

systems = 
[
    ("Пропорциональное H1", [K], [1]),
    ("Интегрирующее H2", [K], [0, T]),
    ("Идеальное дифф. H3", [0, T], [1]),
    ("Реальное дифф. H4", [0, 1], [1, T]),
    ("Первого порядка H5", [K], [1, T]),
]

for a in (0.0, 0.5, 1.5)
    # K / (T^2*s^2 + a*T*s + 1)
    push!(systems, ("Второго порядка H6, a = $a", [K], [1, a * T, T^2]))
end

results = []
for (label, numerator, denominator) in systems
    zero_roots = polynomial_roots(numerator)
    pole_roots = polynomial_roots(denominator)
    push!(results, (label, zero_roots, pole_roots))

    println("$label:")
    println("Нули: ", isempty(zero_roots) ? "нет" : zero_roots)
    println("Полюса: ", isempty(pole_roots) ? "нет" : pole_roots)
end

fig = plot_zeroes_poles(results)

save("poles_zeros_continuous.png", fig)
display(fig)