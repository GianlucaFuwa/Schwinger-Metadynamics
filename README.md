# Schwinger-Metadynamics
Variations of Metadynamics in Simulations of the Schwinger Model (1+1D QED)

## Quick Start
1. Choose parameters of your choice to build a Metapotential first (use: parameters_build.jl to see structure) and do:
```
julia src/run.jl template/parameters_build.jl
```

2. Use built Metapotential to measure observables in second run (use: parameters_build.jl to see structure) and do:
```
julia -t 2 src/jl template/parameters_meas.jl
```
**! Make sure to use 2 or more CPU-Threads when using Parallel Tempering!** 

**! Make sure to use the same Physics and Metapotential parameters for build- and measurement-run !**

Measurements are outputted as .txt files in the chosen directories; you can use them to make plots as you wish.

This code is partly inspired by an old version of the [LatticeQCD.jl](https://github.com/akio-tomiya/LatticeQCD.jl) package by Akio Tomiya et al.
