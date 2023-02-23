# Stochastic population dynamics

This repo contains simulation scripts for modelling the spread of the
covid variants omicron and delta in society, modelled as a Markov jump
process. Three different situations are modelled: societal spread without any
measures implemented, societal spread with quarantine measures, and societal
spread in a vaccinated population.

## Running the scripts

To run the simulations clone the repo and start a Julia session in the base of
the repo. Then run

```julia
include("<script>.jl")
```
where the script is either of `baseline.jl`, `quarantine.jl` or `vaccine.jl`. 
Then run the simulation for a given starting state by entering
```julia
sims = runsim(population_size, I1, I2, endtime, number_of_sims)
```
If you allow Julia several threads, the function will run several simulations in parallel. 

To plot the mean and standard deviation of the simulation, run 
```julia
bandplot(sims, number_of_samples)
```
and to plot the evolution of all simulated states run
```julia
plotter(sims, number_of_samples)
```

