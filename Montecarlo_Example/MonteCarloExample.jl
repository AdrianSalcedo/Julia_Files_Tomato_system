#addprocs()
using DifferentialEquations
using Plots

# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
prob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))

function prob_func(prob,i,repeat)
  ODEProblem(prob.f,rand()*prob.u0,prob.tspan)
end

monte_prob = MonteCarloProblem(prob,prob_func=prob_func)
sim = solve(monte_prob,Tsit5(),num_monte=100)


plotly()
plot(sim,linealpha=0.4)
