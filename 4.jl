using CairoMakie

const SIGNAL_FREQUENCY = 2.0
const PHASE = pi / 6
const OBSERVATION_TIME = 8.0
const TIME_STEP = 0.002
const FREQUENCY_LIMIT = 5.0
const FREQUENCY_STEP = 0.01

times = collect(-OBSERVATION_TIME / 2:TIME_STEP:OBSERVATION_TIME / 2)
frequencies = collect(-FREQUENCY_LIMIT:FREQUENCY_STEP:FREQUENCY_LIMIT)

function cft(signal, times, frequencies)
	weights = fill(times[2] - times[1], length(times))
	weights[[1, end]] ./= 2
	transform = ComplexF64[]

	for frequency in frequencies
		kernel = exp.(-2im * pi * frequency .* times)
		push!(transform, sum(signal .* kernel .* weights))
	end

	return transform
end

function plot_signal_spectrum(title, signal, file_name)
	spectrum = cft(signal, times, frequencies)
	amplitudes = abs.(spectrum)
	phases = angle.(spectrum)

	figure = Figure(size=(1400, 1200))
	axis3 = Axis3(figure[1:2, 1],
		xlabel="Частота, Гц",
		ylabel="Re(X(ν))",
		zlabel="Im(X(ν))",
		title="$title: комплексный спектр",
	)
	lines!(axis3, frequencies, real.(spectrum), imag.(spectrum), linewidth=2)

	amplitude_axis = Axis(figure[3, 1],
		xlabel="Частота, Гц",
		ylabel="|X(ν)|",
		title="$title: амплитудный спектр",
	)
	lines!(amplitude_axis, frequencies, amplitudes, linewidth=2)

	phase_axis = Axis(figure[4, 1],
		xlabel="Частота, Гц",
		ylabel="arg(X(ν)), рад",
		title="$title: фазовый спектр",
	)
	lines!(phase_axis, frequencies, phases, linewidth=2)

	save(file_name, figure)
	display(figure)
end

signals = [
	("x₁(t) = sin(2πft + φ)", sin.(2pi * SIGNAL_FREQUENCY .* times .+ PHASE), "spectrum_x1.png"),
	("x₂(t) = cos(2πft + φ)", cos.(2pi * SIGNAL_FREQUENCY .* times .+ PHASE), "spectrum_x2.png"),
	("x₃(t) = exp(i(2πft + φ))", exp.(1im .* (2pi * SIGNAL_FREQUENCY .* times .+ PHASE)), "spectrum_x3.png"),
]

for (title, signal, file_name) in signals
	plot_signal_spectrum(title, signal, file_name)
end
