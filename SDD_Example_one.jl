#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.SDEProblemLibrary
# load problems
SDEProblemLibrary.importsdeproblems()
prob = SDEProblemLibrary.prob_sde_linear
sol = solve(prob)
