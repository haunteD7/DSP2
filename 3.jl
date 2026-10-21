using Polynomials
using CairoMakie

const SAMPLE_PERIOD = 0.1
const SAMPLE_COUNT = 500

function apply_recurrence(input, numerator, denominator)
    output = zeros(length(input))

    for n in eachindex(input)
        value = 0.0

        for k in eachindex(numerator)
            input_index = n - k + 1
            if input_index >= 1
                value += numerator[k] * input[input_index]
            end
        end

        for k in 2:length(denominator)
            output_index = n - k + 1
            if output_index >= 1
                value -= denominator[k] * output[output_index]
            end
        end

        output[n] = value / denominator[1]
    end

    return output
end

poles_cont = roots(Polynomial([1, 3.2361, 5.2361, 5.2361, 3.2361, 1]))

println("Полюса (непрерывная):")
for pole in poles_cont
    println(pole)
end

b_a = [1.0, 1.0]
a_a = [21.0, -19.0]
b_b = [1.0, 2.0, 1.0]
a_b = [413.361, -798.0, 388.639]
b_c = [1.0, 2.0, 1.0]
a_c = [433.361, -798.0, 368.639]

println("Полюса (дискретная):")
for pole in roots(Polynomial(reverse(a_a)))
    println(pole, " (|z|: ", abs(pole), ")")
end
for pole in roots(Polynomial(reverse(a_b)))
    println(pole, " (|z|: ", abs(pole), ")")
end
for pole in roots(Polynomial(reverse(a_c)))
    println(pole, " (|z|: ", abs(pole), ")")
end


times = collect(0:(SAMPLE_COUNT - 1)) .* SAMPLE_PERIOD
impulse_input = [1.0; zeros(SAMPLE_COUNT - 1)]
response_a = apply_recurrence(impulse_input, b_a, a_a)
response_b = apply_recurrence(impulse_input, b_b, a_b)
response_c = apply_recurrence(impulse_input, b_c, a_c)

impulse_fig = Figure(size=(900, 600))
impulse_ax = Axis(impulse_fig[1, 1],
    xlabel="Время, с",
    ylabel="Амплитуда",
    title="Импульсные характеристики звеньев",
)
lines!(impulse_ax, times, response_a, label="Звено a")
lines!(impulse_ax, times, response_b, label="Звено b")
lines!(impulse_ax, times, response_c, label="Звено c")
axislegend(impulse_ax, position=:rt)
display(impulse_fig)
save("impulse_responses_sections.png", impulse_fig)

overall_b = [
    2.6582840247911e-7,
    1.32914201239555e-6,
    2.6582840247911e-6,
    2.6582840247911e-6,
    1.32914201239555e-6,
    2.6582840247911e-7,
]
overall_a = [
    1.0,
    -4.67669873993035,
    8.75844093261023,
    -8.21010441013895,
    3.85197752757467,
    -0.723606803606714,
]
overall_response = apply_recurrence(impulse_input, overall_b, overall_a)
cascade_a = apply_recurrence(impulse_input, b_a, a_a)
cascade_b = apply_recurrence(cascade_a, b_b, a_b)
composed_response = apply_recurrence(cascade_b, b_c, a_c)

system_fig = Figure(size=(900, 600))
system_ax = Axis(system_fig[1, 1],
    xlabel="Время, с",
    ylabel="Амплитуда",
    title="Импульсная характеристика системы",
)
lines!(system_ax, times, overall_response, label="Заданная система (рекуррентное уравнение)")
lines!(system_ax, times, composed_response, linestyle=:dash,
    label="Композиция звеньев (рекуррентные уравнения)")
axislegend(system_ax, position=:rt)
display(system_fig)
save("impulse_response_system.png", system_fig)
