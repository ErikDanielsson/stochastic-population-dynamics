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

function uniformsample(sim::Vector{Simulation{Int64}}, nsamples, tf)
    nstates = length(sim.states[:, 1])
    tgrid = zeros(nsamples)
    states = zeros(nstates, nsamples)
    tdist = tf / (nsamples - 1)
    tgrid[1] = 0
    states[:, 1] = sim.states[:, 1]
    i = 1
    for j in 2:nsamples
        while sim.t[i+1] < t
            i += 1
        end
    end

end

function uniformsample(sims::Vector{Simulation{Int64}}, nsamples)
    nsims = length(sims)
    nstates = length(sims[1].states[:, 1])
    tgrid = zeros(nsamples)
    states = zeros(nsims, nstates, nsamples)
    tmax = maximum([sim.t[end] for sim in sims])
    tdist = tmax / (nsamples - 1)
    tgrid[1] = 0
    for i in eachindex(sims)
        states[i, :, 1] = sims[i].states[:, 1]
    end
    ind = convert.(Int64, ones(nsims))
    simsinds = Vector(eachindex(sims))
    t = tdist
    for i = 2:nsamples
        maxt = typemin(Float64)
        for j in simsinds
            while sims[j].t[ind[j]+1] < t
                ind[j] += 1
                if sims[j].l == ind[j]
                    popat!(simsinds, findfirst((n) -> n == j, simsinds))
                    break
                end
            end
            if sims[j].t[ind[j]] > maxt
                maxt = sims[j].t[ind[j]]
            end
        end
        tgrid[i] = maxt
        for j in eachindex(sims)
            states[j, :, i] = sims[j].states[:, ind[j]]
        end
        t += tdist
    end
    return tgrid, states
end
