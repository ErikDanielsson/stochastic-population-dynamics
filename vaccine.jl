#=
Subpopulations:
1 S: susceptible population
2 I₁: Infected by virus 1
3 I₂: Infected by virus 2
4 NC₁: Non contagious virus 1
5 NC₂: NOn contagious virus 2
6 M₁: Immune to virus 1
7 M₂: Immune to virus 2

Events: 
Name    Type        Result              Probability/time
1       S -> I₁:    S - 1, I₁ + 1       a₁ * S * I₁ / N
2       S -> I₂:    S - 1, I₂ + 1       a₂ * S * I₂ / N
3       I₁ -> D:    I₁ - 1              b₁ * I₁
4       I₂ -> D:    I₂ - 1              b₂ * I₂
5       I₁ -> NC₁:  I₁ - 1, NC₁ + 1     c₁ * I₁
6       I₂ -> NC₂:  I₂ - 1, NC₂ + 1     c₂ * I₂
7       NC₁ -> D:   NC₁ - 1             n₁ * NC₁
8       NC2 -> D:   NC2 - 1             n₂ * NC₂
9       NC₁ -> M1:  NC₁ - 1, M1 + 1     k₁ * NC₁
10      NC₂ -> M2:  NC₂ - 1, M2 + 1     k2 * NC₂
11      M₁ -> I₂:   M₁ - 1, I₂ + 1      e₁ * M₁ * I₂ / N
12      M₂ -> I₁:   M₂ - 1, I₁ + 1      e₂ * M₂ * I₁ / N 
13      M₁ -> S:    M₁ - 1, S + 1       d₁ * M₁
14      M₂ -> S:    M₂ - 1, S + 1       d₂ * M₂
=#

include("FellerKendall.jl")
include("CustomPlotter.jl")
using LaTeXStrings, ProgressBars

# Define the model graph
const nstates = 7
const graph = [
    (1, 2), # 1
    (1, 3), # 2
    (2, -1), # 3 
    (3, -1), # 4
    (2, 4), # 5 
    (3, 5), # 6
    (4, -1), # 7
    (5, -1), # 8
    (4, 6), # 9
    (5, 7), # 10
    (6, 3), # 11
    (7, 2), # 12
    (6, 1), # 13
    (7, 1), # 14
]

# Define how the transition probabilities depend on the state
getprobs(β, X, N) = @SVector [
    β[1] * X[1] * X[2] / N, # a₁ * S * I₁ / N
    β[2] * X[1] * X[3] / N, # a₂ * S * I₂ / N
    β[3] * X[2], # b₁ * I₁ 
    β[4] * X[3], # b₂ * I₂
    β[5] * X[2], # c₁ * I₁
    β[6] * X[3], # c₂ * I₂
    β[7] * X[4], # n₁ * NC₁
    β[8] * X[5], # n₂ * NC₂
    β[9] * X[4], # k₁ * NC₁
    β[10] * X[5], # k₂ * NC₂
    β[11] * X[6] * X[3] / N, # e₁ * M₁ * I₂ / N
    β[12] * X[7] * X[2] / N, # e₂ * M₂ * I₁ / N
    β[13] * X[6], # d₁ * M₁
    β[14] * X[7], # d₂ * M₂
]


# Real parameters

const totinter = 16 # Total number of interactions per day
const phouse = 2.16 # People per household
const sar_out_d = 0.013 # SAR outside for delta
const sar_out_o = 0.022 # SAR outside for omicron
const sar_in_d = 0.19 # SAR inside for delta
const sar_in_o = 0.30 # SAR inside for omicron
const sickcareful = 0.25 # Corresponds to the carefullness of sick people

const cont_day_d = 6.23 # Number of contagious days for delta
const cont_day_o = 5.35 # Number of contagious days for omicron
const sick_day_d = 10.9 # Number of sick days for delta
const sick_day_o = 9.87 # Number of sick days for omicron
const deadly_d = 0.007 # Baseline deadliness of delta
const deadly_o = 0.004 # Baseline deadliness of omicron
const death_effective_d = 0.9 # Deaths prevented by vaccination for delta
const death_effective_o = 0.59 # Deaths prevented by vaccination for omicron

const meananti_d = 12 * 30 # Average length of antibody protection
const meananti_o = 4 * 30 # Average length of antibody protection

const cross_protection_d = 0.3 # Reduction in the likelihood of omicron infection with delta antibodies
const cross_protection_o = 0.0 # Reduction in the likelihood of delta infection with omicron antibodies

# Transition parameters

const a1 = (totinter - phouse + 1) * sickcareful * sar_out_d + (phouse - 1) * sar_in_d
const a2 = (totinter - phouse + 1) * sickcareful * sar_out_o + (phouse - 1) * sar_in_o

const b1 = (1 - death_effective_d) * deadly_d / sick_day_d
const b2 = (1 - death_effective_o) * deadly_o / sick_day_o

const c1 = (1 - deadly_d) / cont_day_d
const c2 = (1 - deadly_o) / cont_day_o

const d1 = 1 / meananti_d
const d2 = 1 / meananti_o

const e1 = (1 - cross_protection_d) * a2
const e2 = (1 - cross_protection_o) * a1

const n1 = b1
const n2 = b2

const k1 = (1 - n1 * (sick_day_d - cont_day_d)) / (sick_day_d - cont_day_d)
const k2 = (1 - n2 * (sick_day_d - cont_day_d)) / (sick_day_d - cont_day_d)

const β = [
    a1, a2,
    b1, b2,
    c1, c2,
    n1, n2,
    k1, k2,
    e1, e2,
    d1, d2,
]

# Define the transition function after the parameters have been defined
W(X::State, N::Int) = getprobs(β, X, N)
const model = Model(nstates, graph, W)

# Now we can define a problem instance
getX0(N, I1, I2) = State{nstates}([N - I1 - I2, I1, I2, 0, 0, 0, 0])


function runsim(N, I1, I2, tf, nsims::Int)
    sims = Vector{Simulation{Int64}}(undef, nsims)
    Threads.@threads for i in ProgressBar(1:nsims)
        sims[i] = simulation(getX0(N, I1, I2), model, tf)
    end
    return sims
end

const labels = [L"S" L"D" L"I_1" L"I_2" L"NC_1" L"NC_2" L"M_1" L"M_2"]
bandplot(sims::Vector, nsamples) = bandplot(sims, nsamples, nstates, labels)
plotter(sim::Simulation, inc) = plotter(sim, inc, nstates, labels)
plotter!(sim::Simulation, inc) = plotter!(sim, inc, nstates, labels)
plotter(sims::Vector{Simulation}, inc) = plotter(sims, inc, nstates, labels)