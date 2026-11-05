using CairoMakie

const T = 0.5
const K = 1.0
const TAU = 2.0
const A_VALUES = (0.0, 0.5, 1.5)
const OMEGA = 10 .^ range(-2, 2, length=1600)

function h1(omega, a=0.0)
    return fill(complex(K), length(omega))
end

function h2(omega, a=0.0)
    return K ./ (im .* T .* omega)
end

function h3(omega, a=0.0)
    return im .* T .* omega
end

function h4(omega, a=0.0)
    return im .* omega ./ (1 .+ im .* T .* omega)
end

function h5(omega, a=0.0)
    return K ./ (1 .+ im .* T .* omega)
end

function h6(omega, a)
    return K ./ (1 .- (T .* omega).^2 .+ im .* a .* T .* omega)
end

function h7(omega, a=0.0)
    return K .* exp.(-im .* TAU .* omega)
end

function h6_cutoff(a)
    if a == 0
        return nothing
    end

    # |H₆(jω)| = K / sqrt(2), r = Tω.
    r_squared = (2 - a^2 + sqrt((a^2 - 2)^2 + 4)) / 2
    return (sqrt(r_squared) / T, 20 * log10(K) - 3)
end

function h4_cutoff(a=0.0)
    return (1 / T, 20 * log10(1 / T) - 3)
end

function h5_cutoff(a=0.0)
    return (1 / T, 20 * log10(K) - 3)
end

systems = [
    ("H₁: пропорциональное звено", h1, (nothing,), nothing, "h1"),
    ("H₂: интегрирующее звено", h2, (nothing,), nothing, "h2"),
    ("H₃: идеальное дифференцирующее звено", h3, (nothing,), nothing, "h3"),
    ("H₄: реальное дифференцирующее звено", h4, (nothing,), h4_cutoff, "h4"),
    ("H₅: звено первого порядка", h5, (nothing,), h5_cutoff, "h5"),
    ("H₆: звено второго порядка", h6, A_VALUES, h6_cutoff, "h6"),
    ("H₇: идеальное запаздывающее звено", h7, (nothing,), nothing, "h7"),
]

function unwrap_phase(values)
    result = copy(angle.(values))
    phase_shift = 0.0

    for index in 2:length(result)
        difference = result[index] - result[index - 1]
        if difference > pi
            phase_shift -= 2pi
        elseif difference < -pi
            phase_shift += 2pi
        end
        result[index] += phase_shift
    end

    return result
end

function logarithmic_amplitude(values)
    return 20 .* log10.(abs.(values))
end

function plot_system(system)
    title, transfer_function, parameters, cutoff_function, filename = system
    amplitude_figure = Figure(size=(950, 600))
    amplitude_axis = Axis(amplitude_figure[1, 1],
        xscale=log10,
        xlabel="ω, рад/с",
        ylabel="20 log₁₀(A), дБ",
        title="$title: ЛАЧХ",
    )

    phase_figure = Figure(size=(950, 600))
    phase_axis = Axis(phase_figure[1, 1],
        xscale=log10,
        xlabel="ω, рад/с",
        ylabel="Φ(ω), рад",
        title="$title: ФЧХ",
    )

    for parameter in parameters
        value = parameter === nothing ? 0.0 : parameter
        transfer = transfer_function(OMEGA, value)
        label = parameter === nothing ? title : "a = $parameter"

        lines!(amplitude_axis, OMEGA, logarithmic_amplitude(transfer),
            linewidth=2, label=label)
        lines!(phase_axis, OMEGA, unwrap_phase(transfer),
            linewidth=2, label=label)

        if cutoff_function !== nothing
            cutoff = cutoff_function(value)
            if cutoff !== nothing
                omega_cutoff, cutoff_db = cutoff
                vlines!(amplitude_axis, [omega_cutoff], linestyle=:dash, linewidth=1)
                scatter!(amplitude_axis, [omega_cutoff], [cutoff_db], markersize=9,
                    label="ωср = $(round(omega_cutoff, digits=3))")
                hlines!(amplitude_axis, [cutoff_db], color=:black, linestyle=:dash,
                    linewidth=1, label="уровень -3 дБ")
            end
        end
    end

    if filename == "h7"
        ylims!(amplitude_axis, -1, 1)
    end

    axislegend(amplitude_axis, position=:lb)
    axislegend(phase_axis, position=:lb)
    save("$(filename)_lachx.png", amplitude_figure)
    save("$(filename)_fchx.png", phase_figure)
end

for system in systems
    plot_system(system)
end
