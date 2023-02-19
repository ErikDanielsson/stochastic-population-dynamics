#=
Subpopulations:
1 S: susceptible population
2 I₁: Infected by virus 1
3 I₂: Infected by virus 2
4 M₁: Immune to virus 1
5 M₂: Immune to virus 2

Events: 
Name    Type        Result              Probability/time
1       S -> E₁:    S - 1, I₁ + 1       a₁ * S * I₁ / N
2       S -> E₂:    S - 1, I₂ + 1       a₂ * S * I₂ / N
3       I₁ -> D:    I₁ - 1              b₁ * I₁
4       I₂ -> D:    I₂ - 1              b₂ * I₂
5       I₁ -> M₁:   I₁ - 1, M₁ + 1      d₁ * I₁
6       I₂ -> M₂:   I₂ - 1, M₂ + 1      d₂ * I₂
7       M₁ -> S:    M₁ - 1, S + 1       f₁ * M₁
8       M₂ -> S:    M₂ - 1, S + 1       f₂ * M₂
9       M₁ -> I₂:   M₁ - 1, I₂ + 1      e₁ * M₁ * I₂ / N
10      M₂ -> I₁:   M₂ - 1, I₁ + 1      e₂ * M₂ * I₁ / N 
=#

include("FellerKendall.jl")
include("CustomPlotter.jl")
using LaTeXStrings, ProgressBars

# Define the model graph
const nstates = 5
const graph = [
    (1, 2),
    (1, 3),
    (2, -1),
    (3, -1),
    (2, 4),
    (3, 5),
    (4, 1),
    (5, 1),
    (4, 3),
    (5, 2),
]

# Define how the transition probabilities depend on the state
getprobs(β, X, N) = @SVector [
    β[1] * X[1] * X[2] / N, # Infection 1
    β[2] * X[1] * X[3] / N, # Infection 2
    β[3] * X[2], # Death 1
    β[4] * X[3], # Death 2
    β[5] * X[2], # Immune 1
    β[6] * X[3], # Immune 2
    β[7] * X[4], # Susceptible 1
    β[8] * X[5], # Susceptible 2
    β[9] * X[4] * X[3] / N, # crossinfection (δ -> o)
    β[10] * X[5] * X[2] / N, # Crossinfectiona (o -> δ) 
]

# Real parameters

const totinter = 16 # Total number of interactions per day
const phouse = 2.16 # People per household
const sar_out_d = 0.046 # SAR outside for delta 
const sar_out_o = 0.050 # SAR outside for omicron
const sar_in_d = 0.27 # SAR inside for delta
const sar_in_o = 0.28 # SAR inside for omicron
const sickcareful = 0.25 # Corresponds to the carefullness of sick people

const cont_day_d = 6.23 # Number of contagious days for delta
const cont_day_o = 5.35 # Number of contagious days for omicron
const sick_day_d = 10.9 # Number of sick days for delta
const sick_day_o = 9.87 # Number of sick days for omicron
const deadly_d = 0.007 # Baseline deadliness of delta
const deadly_o = 0.004 # Baseline deadliness of omicron

const meananti_d = 12 * 30 # Average length of antibody protection
const meananti_o = 4 * 30 # Average length of antibody protection

const cross_protection_d = 0.3 # Reduction in the likelihood of omicron infection with delta antibodies
const cross_protection_o = 0.277 # Reduction in the likelihood of delta infection with omicron antibodies

# Transition parameters

const a1 = (totinter - phouse + 1) * sickcareful * sar_out_d + (phouse - 1) * sar_in_d
const a2 = (totinter - phouse + 1) * sickcareful * sar_out_o + (phouse - 1) * sar_in_o

const b1 = deadly_d / sick_day_d
const b2 = deadly_o / sick_day_o

const c1 = (1 - deadly_d) / cont_day_d
const c2 = (1 - deadly_o) / cont_day_o

const d1 = 1 / meananti_d
const d2 = 1 / meananti_o

const e1 = (1 - cross_protection_d) * a1
const e2 = (1 - cross_protection_o) * a2

const β = [
    a1, a2,
    b1, b2,
    c1, c2,
    d1, d2,
    e1, e2,
]

# Define the transition function after the parameters have been defined
W(X::State, N::Int) = getprobs(β, X, N)
const model = Model(nstates, graph, W)

# Now we can define a problem instance
getX0(N, I1, I2) = State{nstates}([N - I1 - I2, I1, I2, 0, 0])

function runsim(N, I1, I2, tf, nsims::Int)
    # Convenience function for running the simulation on several threads
    sims = Vector{Simulation{Int64}}(undef, nsims)
    Threads.@threads for i in Progressbar(1:nsims)
        sims[i] = simulation(getX0(N, I1, I2), model, tf)
    end
    return sims
end

# Define the plotting functions for the problem instance
const labels = [L"S" L"D" L"I_1" L"I_2" L"M_1" L"M_2"]
bandplot(sims::Vector, nsamples) = bandplot(sims, nsamples, nstates, labels)
plotter(sim::Simulation, inc) = plotter(sim, inc, nstates, labels)
plotter!(sim::Simulation, inc) = plotter!(sim, inc, nstates, labels)
plotter(sims::Vector{Simulation{Int}}, inc) = plotter(sims, inc, nstates, labels)