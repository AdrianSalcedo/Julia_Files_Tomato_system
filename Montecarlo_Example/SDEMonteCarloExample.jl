using DifferentialEquations
using Plots
#using Plotly

function f(du,u,p,t)
  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
end

function g(du,u,p,t)
  du[1] = p[3]*u[1]
  du[2] = p[4]*u[2]
end

p = [1.5,1.0,0.1,0.1]
prob = SDEProblem(f,g,[1.0,1.0],(0.0,10.0),p)

function prob_func(prob,i,repeat)
  prob.p[3:4] = 0.3rand(2)
  prob
end

monte_prob = MonteCarloProblem(prob,prob_func=prob_func)
sim = solve(monte_prob,SRIW1(),trajectories=10)
plotly()
plot(sim,linealpha=0.6,color=:blue,vars=(0,1),title="Phase Space Plot")
plot!(sim,linealpha=0.6,color=:red,vars=(0,2),title="Phase Space Plot")

summ = MonteCarloSummary(sim,0:0.1:10)
pyplot() # Note that plotly does not support ribbon plots
plot(summ,fillalpha=0.5)
#A=sim.u acceso a todas las simulaciones. Cada simulacion contiene  u , t, NoiseProcces (W)
#A[1]  es la primera simulacion del systema.
