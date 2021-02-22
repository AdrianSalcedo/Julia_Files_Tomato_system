using DifferentialEquations
α=1
β=1
u₀=1/2
f(u,p,t) = α*u
g(u,p,t) = β*u
dt = 1//2^(4)
tspan = (0.0,1.0)

prob = SDEProblem(f,g,u₀,(0.0,1.0))
sol = solve(prob,SRIW1(),dt=dt,adaptive=false)
using Plots; plotly() # Using the Plotly backendPkg.status("DifferentialEquations")
plot(sol)

ensembleprob = EnsembleProblem(prob)

sol = solve(ensembleprob,EnsembleThreads(),trajectories=1000)

using DifferentialEquations.EnsembleAnalysis
summ = EnsembleSummary(sol,0:0.01:1)
plot(summ,labels="Middle 95%")
summ = EnsembleSummary(sol,0:0.01:1;quantiles=[0.25,0.75])
plot!(summ,labels="Middle 50%",legend=true)
