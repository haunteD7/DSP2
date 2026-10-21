using Polynomials

include("plot.jl")

m = 100
K = 1

function polynomial_roots(coefficients)
    return roots(Polynomial(coefficients))
end

systems = 
[
    ("Пропорциональное H1", [K], [1]),
    ("Интегрирующее H2", [0, K], [-m, m]),
    ("Дифф. H3", [-K*m, K*m], [0, 1]),
    ("Первого порядка H5", [K], [m + 1, 1]),
    ("Второго порядка H6", [K], [1, m, m^2]),
    ("Идеальное запаздывающее H7", [K], [zeros(m); 1]),
    
]

results = []
for (label, num_coeffs, den_coeffs) in systems
    zero_roots = polynomial_roots(num_coeffs)
    pole_roots = polynomial_roots(den_coeffs)

    push!(results, (label, zero_roots, pole_roots))
end


fig = plot_zeroes_poles(results)
save("poles_zeros_discrete.png", fig)