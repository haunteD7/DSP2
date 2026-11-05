using CairoMakie
using FFTW

const SAMPLE_RATE = 200.0
const INITIAL_DURATION = 2.0
const EXTENDED_DURATION = 8.0
const INITIAL_FREQUENCIES = (1.0, 10.0, 30.0)
const SHIFTED_FREQUENCIES = (1.25, 10.0, 30.0)
const AMPLITUDES = (1.0, 0.5, 0.3)
const ZERO_PADDED_SIZE = 4096
const EXTENDED_FFT_SIZE = 16384

function create_signal(frequencies, amplitudes, duration)
    sample_count = round(Int, SAMPLE_RATE * duration)
    times = (0:sample_count - 1) ./ SAMPLE_RATE
    values = zeros(sample_count)

    for (frequency, amplitude) in zip(frequencies, amplitudes)
        values .+= amplitude .* sin.(2pi * frequency .* times)
    end

    return values
end

function zero_pad(values, target_length)
    return [values; zeros(target_length - length(values))]
end

function calculate_spectrum(values, fft_length; normalization_length=length(values))
    padded_values = zero_pad(values, fft_length)
    transform = fft(padded_values)
    frequencies = (0:fft_length ÷ 2) .* SAMPLE_RATE ./ fft_length
    amplitudes = 2 .* abs.(transform[1:fft_length ÷ 2 + 1]) ./ normalization_length
    amplitudes[1] /= 2

    return frequencies, amplitudes
end

function plot_spectrum(frequencies, amplitudes, title, filename; points=false)
    figure = Figure(size=(950, 600))
    axis = Axis(figure[1, 1],
        xlabel="f, Гц",
        ylabel="Амплитуда",
        title=title,
        limits=(0, 50, nothing, nothing),
    )
    if points
        scatter!(axis, frequencies, amplitudes, markersize=7)
    else
        lines!(axis, frequencies, amplitudes, linewidth=2)
    end
    save(filename, figure)
end

function hann_window(window_length)
    indices = 0:window_length - 1
    return 0.5 .* (1 .- cos.(2pi .* indices ./ (window_length - 1)))
end

function blackman_window(window_length)
    indices = 0:window_length - 1
    return 0.42 .- 0.5 .* cos.(2pi .* indices ./ (window_length - 1)) .+
        0.08 .* cos.(4pi .* indices ./ (window_length - 1))
end

function I_0(value)
    term = 1.0
    result = 1.0

    for index in 1:100
        term *= (value^2 / 4) / index^2
        result += term
        if abs(term) < eps(result)
            break
        end
    end

    return result
end

function kaiser_window(window_length, beta=5.0)
    indices = 0:window_length - 1
    argument = beta .* sqrt.(1 .- (2 .* indices ./ (window_length - 1) .- 1).^2)
    return I_0.(argument) ./ I_0(beta)
end

experiments = [
    ("Исходный сигнал: f = 1, 10, 30 Гц", INITIAL_FREQUENCIES,
        INITIAL_DURATION, 400, "7_initial.png", true),
    ("Первая компонента смещена на 1.25 Гц", SHIFTED_FREQUENCIES,
        INITIAL_DURATION, 400, "7_shifted.png", true),
    ("Дополнение нулями до N = 4096", INITIAL_FREQUENCIES,
        INITIAL_DURATION, ZERO_PADDED_SIZE, "7_zero_padding.png", false),
    ("Сигнал длительностью 8 с, N = 16384", INITIAL_FREQUENCIES,
        EXTENDED_DURATION, EXTENDED_FFT_SIZE, "7_extended.png", false),
]

for (title, frequencies, duration, fft_length, filename, points) in experiments
    values = create_signal(frequencies, AMPLITUDES, duration)
    spectrum_frequencies, spectrum_amplitudes = calculate_spectrum(values, fft_length)
    plot_spectrum(spectrum_frequencies, spectrum_amplitudes, title, filename, points=points)
end

windows = [
    ("Без окна", window_length -> ones(window_length)),
    ("Ханна", hann_window),
    ("Кайзера, β = 5", kaiser_window),
    ("Блэкмана", blackman_window),
]

windowed_values = create_signal(INITIAL_FREQUENCIES, AMPLITUDES, INITIAL_DURATION)
window_figure = Figure(size=(1100, 850))

for (index, (title, window_function)) in enumerate(windows)
    window = window_function(length(windowed_values))
    spectrum_frequencies, spectrum_amplitudes = calculate_spectrum(
        windowed_values .* window,
        ZERO_PADDED_SIZE,
        normalization_length=sum(window),
    )
    axis = Axis(window_figure[(index - 1) ÷ 2 + 1, (index - 1) % 2 + 1],
        xlabel="f, Гц",
        ylabel="Амплитуда",
        title=title,
        limits=(0, 50, nothing, nothing),
    )
    lines!(axis, spectrum_frequencies, spectrum_amplitudes, linewidth=2)
end

save("7_windows.png", window_figure)
