#=
The Feller-Kendall algorithm for a Markov jump processes

Author: Erik Danielsson

=#
using Distributions, StaticArrays

struct Model
    nstates::Int64
    graph::Vector{Tuple{Int64,Int64}}
    W
end

"""
Compile the model by recasting the graph tuples as their corresponding
actions on the state vector
"""
function compilemodel(model::Model)::OptModel
    nevents = length(model.graph)
    nstates = model.nstates
    actions = zeros(Int64, nstates, nevents)
    for (i, edge) in enumerate(model.graph)
        action = zeros(nstates)
        s, e = edge
        if s != -1
            action[s] = -1
        end
        if e != -1
            action[e] = 1
        end
        actions[:, i] = action
    end
    return OptModel{nstates,nevents}(
        SMatrix{nstates,nevents,Int64}(actions),
        model.W
    )
end

struct OptModel{nstates,nevents}
    actions::SMatrix{nstates,nevents,Int64}
    W
end

State{N} = SVector{N,Int64}

struct Simulation{T}
    l::Int64
    t::Vector{Float64}
    states::Matrix{T}
    popsize::T
end

function Simulation{Int64}(tgrid::Vector{Float64},
    states::Vector{State},
    popsize::Int64
)
    states = reduce(hcat, states) # Matrix conversion
    deaths = popsize .- sum(states, dims=1) # Calculate death count

    mstates = [states[1, :]'; deaths; states[2:end, :]]# vcat(states, deaths)
    Simulation(length(tgrid), tgrid, mstates, popsize)
end

"""
This is the main simulation step: calculating transition probabilities and
sampling the distributions. 
"""
function nextevent(X::State, popsize::Int64, optmodel::OptModel)
    # Compute the event probabilities for the current state
    W = optmodel.W(X, popsize)
    # Sample the next event time
    R = sum(W)
    # If all probabilities are 0, the we should end the simulation, 
    # so we send an exit signal
    if R == 0.0
        return -1, X
    end
    # Sample the next time from an exponential distribution with mean R
    Δt = rand(Exponential(1 / R))
    # Sample a uniform distribution to get the next event
    q = rand(Uniform(0, R))
    for i in eachindex(W)
        # Choose an event based on the sizes of the probabilities 
        # The effect is determined by the compiled problem graph
        if q < W[i]
            nextX = X + optmodel.actions[:, i]
            return Δt, nextX
        end
        q -= W[i]
    end
end

"""
The Feller-Kendall algorithm
"""
function simulation(X0::State, model::Model, tf::Number)
    optmodel = compilemodel(model)
    N = sum(X0)
    X = X0
    t = 0
    states = Vector{State}()
    tgrid = Vector{Float64}()
    push!(states, X)
    push!(tgrid, 0)
    while t < tf
        Δt, X = nextevent(X, N, optmodel)
        if Δt == -1 # Exit signal
            break
        end
        t += Δt
        push!(states, X)
        push!(tgrid, t)
    end
    return Simulation{Int64}(tgrid, states, N)
end

"""
If we have several simulations, we need to join them. This means defining each 
simulation vector for each time point. 

The extra parameter is for only saving every nth point. Useful for plotting
"""
function joinsims(sims::Vector{Simulation{Int64}}, inc::Int64=1)
    nstates = length(sims[1].states[:, 1]) # Number of states in the current problem
    nsims = length(sims)
    ntimes = sum([length(s.t) for s in sims]) - (nsims - 1) # All simulations are equal for t = 0, hence the subtraction
    nsteps = div(ntimes, inc) + 1
    tgrid = zeros(nsteps)
    states = zeros(nsims, nstates, nsteps)
    ind = convert.(Int64, ones(nsims)) # Current indices for each simulation
    # Next we need to create a vector containing indices for the simulations
    # we are currently summing. This is to ensure that we don't need to 
    # compare a simulation that ended earlier
    simsinds = Vector(eachindex(sims))
    tgrid[1] = 0
    for (i, sim) in enumerate(sims)
        states[i, :, 1] = sims[i].states[:, 1]
    end
    for k in 2:ntimes
        j = -1
        nextt = typemax(Float64)
        for i in simsinds[1:end]
            if sims[i].l ≤ ind[i]
                popat!(simsinds, findfirst((n) -> n == i, simsinds))
                continue
            end
            if sims[i].t[ind[i]+1] < nextt
                nextt = sims[i].t[ind[i]+1]
                j = i
            end
        end
        ind[j] += 1
        if k % inc == 1
            tgrid[div(k, inc)+1] = nextt
            for (i, (sim, j)) in enumerate(zip(sims, ind))
                states[i, :, div(k, inc)+1] = sim.states[:, j]
            end
        end
    end
    return tgrid, states # Simulation{Float64}(length(tgrid), tgrid, states, sims[1].popsize)
end
