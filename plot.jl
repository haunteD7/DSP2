using CairoMakie

function plot_zeroes_poles(results)
    fig = Figure(size=(900, 650))
    ax = Axis(fig[1, 1],
        xlabel="Re(s)",
        ylabel="Im(s)",
        title="Нули и полюса систем",
    )
    vlines!(ax, [0], color=:gray, linewidth=1)
    hlines!(ax, [0], color=:gray, linewidth=1)


    for (label, zero_roots, pole_roots) in results
        has_label = false

        if !isempty(zero_roots)
            scatter!(ax, real.(zero_roots), imag.(zero_roots),
                marker=:circle, markersize=14,
                label=label,
            )
            has_label = true
        end

        if !isempty(pole_roots)
            scatter!(ax, real.(pole_roots), imag.(pole_roots),
                marker=:xcross, markersize=14,
                label=has_label ? nothing : label,
            )
        end
    end

    axislegend(ax, position=:lt)

    return fig
end