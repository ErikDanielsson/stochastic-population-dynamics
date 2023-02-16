using Distributions, Plots

function baseplotter(plotfunc, sim::Simulation, inc::Int, nstates, labels)
    plotfunc(sim.t[1:inc:end],
        sim.states'[1:inc:end, :],
        layout=(ceil(Int64, (nstates + 1) / 2), 2),
        title=labels,
        legend=false,
    )
end

function plotter(sim::Simulation, inc::Int, nstates, labels)
    baseplotter(plot, sim, inc, nstates, labels)
end

function plotter!(sim::Simulation, inc::Int, nstates, labels)
    baseplotter(plot!, sim, inc, nstates, labels)
end

function plotter(sims::Vector{Simulation{Int}}, inc::Int, nstates, labels)
    plotter(sims[1], inc, nstates, labels)
    for s in sims[2:end]
        plotter!(s, inc, nstates, labels)
    end
    current()
end

function bandplot(sims::Vector, nsamples, nstates, labels)
    tgrid, states = uniformsample(sims, nsamples)
    meanstates = mean(states, dims=1)
    σ = sqrt.(dropdims(var(states, mean=meanstates, dims=1), dims=1))
    μ = dropdims(meanstates, dims=1)
    min = dropdims(minimum(states, dims=1), dims=1)
    max = dropdims(maximum(states, dims=1), dims=1)

    layout = (ceil(Int64, (nstates + 1) / 2), 2)
    plot(tgrid,
        min',
        fillrange=max',
        fillalpha=0.15,
        linealpha=0,
        c=:blue,
        title=labels,
        layout=layout,
        legend=false
    )


    plot!(tgrid,
        (μ + σ)',
        fillrange=(μ - σ)',
        fillalpha=0.35,
        linealpha=0,
        c=:blue,
        title=labels,
        layout=layout,
        legend=false
    )

    plot!(tgrid,
        μ',
        c=:blue,
        title=labels,
        layout=layout,
        legend=false,
    )
end

function uniformsample(sim::Simulation{Int64}, tdist, nsamples, nstates)
    tgrid = zeros(nsamples)
    sampstates = zeros(nstates, nsamples)
    tgrid[1] = 0
    simstates = load_sim_states(sim)
    t = tdist
    j = 1
    for i = 2:nsamples
        while sim.t[j+1] < t
            j += 1
            if sim.l == j
                break
            end
        end
        tgrid[i] = sim.t[j]
        if sim.l == j
            sampstates[:, i:end] .= simstates[:, j]
            break
        else
            sampstates[:, i] = simstates[:, j]
        end
        t += tdist
    end
    return tgrid, sampstates
end

function uniformsample(sims::Vector{Simulation{Int64}}, nsamples)
    nsims = length(sims)
    nstates = sims[1].nstates
    tgrid = zeros(nsamples)
    states = zeros(nsims, nstates, nsamples)
    tmax = maximum([sim.t[end] for sim in sims])
    tdist = tmax / (nsamples - 1)
    for i in 1:nsims
        sampgrid, sampstates = uniformsample(sims[i], tdist, nsamples, nstates)
        states[i, :, :] = sampstates
        println("a ", tgrid)
        println("b ", sampgrid)
        tgrid = max.(tgrid, sampgrid)
        println("c ", tgrid)
    end
    return tgrid, states
end

