using CairoMakie

const M = 100
const K = 1.0
const SAMPLE_RATES = (100.0, 10000.0)
const OMEGA_MIN = 1e-4
const OMEGA = vcat(
    collect(range(-pi, -OMEGA_MIN, length=1000)),
    collect(range(OMEGA_MIN, pi, length=1000)),
)
const INPUT_SIGNALS = [
    ("x₁(t) = sin(2π · 0.1t)", 0.1, 100.0),
    ("x₂(t) = cos(2π · 1t)", 1.0, 100.0),
    ("x₃(t) = exp(i(2π · 10t))", 10.0, 10000.0),
]

function h1(omega)
    z = exp.(im .* omega)
    return fill(complex(K), length(z))
end

function h2(omega)
    z = exp.(im .* omega)
    return K .* z ./ (M .* (z .- 1))
end

function h3(omega)
    z = exp.(im .* omega)
    return K .* M .* (z .- 1) ./ z
end

function h5(omega)
    z = exp.(im .* omega)
    return K .* M ./ (z .+ (M + 1))
end

function h6(omega)
    z = exp.(im .* omega)
    return K ./ (M^2 .* z.^2 .+ M .* z .+ 1)
end

function h7(omega)
    z = exp.(im .* omega)
    return K .* z .^ (-M)
end

function logarithmic_amplitude(values)
    return 20 .* log10.(abs.(values))
end

function plot_system(system)
    title, transfer_function, stable, filename = system
    amplitude_figure = Figure(size=(1400, 850))
    phase_figure = Figure(size=(1400, 850))

    amplitude_axes = Axis[]
    phase_axes = Axis[]
    axis_titles = [
        "ω, рад/отсчёт",
        "f, Гц (fₛ = 100 Гц)",
        "f, Гц (fₛ = 10 кГц)",
    ]

    for column in 1:3
        amplitude_axis = Axis(amplitude_figure[1, column],
            xlabel=axis_titles[column],
            ylabel=column == 1 ? "20 log₁₀(A), дБ" : "",
            title=column == 1 ? "$title: ЛАЧХ" : "",
        )
        phase_axis = Axis(phase_figure[1, column],
            xlabel=axis_titles[column],
            ylabel=column == 1 ? "Φ(ω), рад" : "",
            title=column == 1 ? "$title: ФЧХ" : "",
        )
        push!(amplitude_axes, amplitude_axis)
        push!(phase_axes, phase_axis)
    end

    transfer = transfer_function(OMEGA)
    phase_omega = title == "H₇" ? collect(range(-3pi / M, 3pi / M, length=1600)) : OMEGA
    phase_transfer = transfer_function(phase_omega)
    amplitudes = logarithmic_amplitude(transfer)
    phases = angle.(phase_transfer)

    for column in 1:3
        frequencies = column == 1 ? OMEGA : OMEGA .* SAMPLE_RATES[column - 1] ./ (2pi)
        phase_frequencies = column == 1 ? phase_omega : phase_omega .* SAMPLE_RATES[column - 1] ./ (2pi)
        lines!(amplitude_axes[column], frequencies, amplitudes,
            linewidth=2, label=title)
        lines!(phase_axes[column], phase_frequencies, phases,
            linewidth=2, label=title)
    end

    for axis in amplitude_axes
        axislegend(axis, position=:lb)
    end
    for axis in phase_axes
        axislegend(axis, position=:lb)
    end

    if title == "H₇"
        for axis in amplitude_axes
            ylims!(axis, -1, 1)
        end
    end

    save("$(filename)_discrete_lachx.png", amplitude_figure)
    save("$(filename)_discrete_fchx.png", phase_figure)
end

function print_output_characteristics()
    for (signal_title, input_frequency, sample_rate) in INPUT_SIGNALS
        omega = 2pi * input_frequency / sample_rate
        println()
        println("$signal_title, fₛ = $sample_rate Гц, ω = $omega рад/отсчёт")

        for (title, transfer_function, stable, filename) in systems
            if !stable
                continue
            end

            transfer = transfer_function([omega])[1]
            output_amplitude = abs(transfer)
            output_phase = angle(transfer)
            println("$title: A = $output_amplitude, Φ = $output_phase рад")
        end
    end
end

systems = [
    ("H₁", h1, true, "h1"),
    ("H₂", h2, false, "h2"),
    ("H₃", h3, true, "h3"),
    ("H₅", h5, false, "h5"),
    ("H₆", h6, true, "h6"),
    ("H₇", h7, true, "h7"),
]

for system in systems
    plot_system(system)
end

print_output_characteristics()
