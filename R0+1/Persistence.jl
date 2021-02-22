using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using IterableTables, DataFrames, DataTables
#using DiffEqGPU
#using CuArrays
using StochasticDiffEq

beta_p = 0.8
r_1 = 0.6
r_2 = 0.6
b = 0.5
beta_v = 0.8
theta = 0.818469652
mu = 0.058596601
gamma = 0.4
sigma_L = 0.378173863
sigma_I = 0.030015876
sigma_v = 0.254146903
N_v = mu/gamma
r= max(r_1,r_2)
u_0 = [97.0,1.0,2.0,3.0,4.0]
T = 100.0
time = (0.0,T)
N_p = u_0[1]+u_0[2]+u_0[3]
dt=0.01
t_s=range(0.0,T, step=1.0)


Rs0 = beta_p*beta_v/(r*gamma)
print("Rs0=",Rs0)

function F_Drift(du,u,p,t)
     @inbounds begin
        du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
        du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
        du[3] = b*u[2]-r_2*u[3]
        du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
        du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
    end
    nothing

end

function G_Diffusion(du,u,p,t)
    @inbounds begin
        du[1,1] = sigma_L*u[2]*u[1]/N_p+sigma_I*u[1]*u[3]/N_p
        du[1,2] = 0
        du[2,1] = -sigma_L*u[1]*u[2]/N_p
        du[2,2] = 0
        du[3,1] = -sigma_I*u[1]*u[3]/N_p
        du[3,2] = 0
        du[4,1] = 0
        du[4,2] = -sigma_v*u[4]
        du[5,1] = 0
        du[5,2] = -sigma_v*u[5]
    end
    nothing
end

########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,noise_rate_prototype=zeros(5,2))
sol = solve(prob_sde_tomato_sys)

#################w###############################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################
monte_prob = MonteCarloProblem(prob_sde_tomato_sys)
sim = solve(monte_prob,trajectories=10)
summ_1 = MonteCarloSummary(sim,t_s)
summ_2 = MonteCarloSummary(sim,t_s;quantiles=[0.25,0.75])
################################################################################
################### To comprovate the conservation law #########################
################################################################################
#=
Datos2=DataFrame()
component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
component = transpose(component) #transpose to obtain any*5 data matrix
component = vcat(component...) #to obtain shape for dataframe
component = vcat(component...) # again do a reshape
variables = DataFrame(component) # define first data frame
Datos1 = DataFrame(S_p = variables[:,1], L_p = variables[:,2], I_p = variables[:,3]) #only some variables
Datos2 = append!(Datos2, Datos1)
=#

################################################################################
#############################    Plot Summary   ################################
################################################################################

plotly()
title = plot(title = "R_s =$Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
q1=plot(summ_1,idxs=(1),labels="Middle 95%")
q1=plot!(summ_2,idxs=(1),labels="Middle 50%",legend=true)
q2=plot(summ_1,idxs=(2),labels="Middle 95%")
q2=plot!(summ_2,idxs=(2),labels="Middle 50%",legend=true)
q3=plot(summ_1,idxs=(3),labels="Middle 95%")
q3=plot!(summ_2,idxs=(3),labels="Middle 50%",legend=true)
q4=plot(summ_1,idxs=(4),labels="Middle 95%")
q4=plot!(summ_2,idxs=(4),labels="Middle 50%",legend=true)
q5=plot(summ_1,idxs=(5),labels="Middle 95%")
q5=plot!(summ_2,idxs=(5),labels="Middle 50%",legend=true)

plot(q1,q2,q3,q4,q5,title,layout = @layout([ [A B C]; [D E F]]), label="")

#=
################################################################################
########################    data  Media    #####################################
################################################################################
time = summ_1.t
xu = summ_1.u
xu_glued = hcat(xu...)
Xu1 = xu_glued[1:5:end]
Xu3 = xu_glued[3:5:end]
Xu5 = xu_glued[5:5:end]

DF1 = DataFrame(t = time,S_p = Xu1,I_p = Xu3, I_v = Xu5)
CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Persistence_Mean.csv",DF1)

################################################################################
######################  quartile default data  Media ###########################
################################################################################
time = summ_1.t
xqlow = summ_1.qlow.u
X_glued = hcat(xqlow...)
X1 = X_glued[1:5:end]
X3 = X_glued[3:5:end]
X5 = X_glued[5:5:end]

yqhigh = summ_1.qhigh.u
Y_glued = hcat(yqhigh...)
Y1 = Y_glued[1:5:end]
Y3 = Y_glued[3:5:end]
Y5 = Y_glued[5:5:end]

DF2 = DataFrame(t = time,S_p = X1, I_p = X3, I_v = X5)
DF3 = DataFrame(t = time,S_p = Y1, I_p = Y3, I_v = Y5)
CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Persistence_Mean_Quartile_low.csv",DF2)
CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Persistence_Mean_Quartile_high.csv",DF3)


################################################################################
#######################   Data quantiles 50%    ################################
################################################################################
time = summ_2.t
xxqlow = summ_2.qlow.u
XX_glued = hcat(xxqlow...)
XX1 = XX_glued[1:5:end]
XX3 = XX_glued[3:5:end]
XX5 = XX_glued[5:5:end]

yyqhigh = summ_2.qhigh.u
YY_glued = hcat(yyqhigh...)
YY1 = YY_glued[1:5:end]
YY3 = YY_glued[3:5:end]
YY5 = YY_glued[5:5:end]

DF4 = DataFrame(t = time,S_p = XX1, I_p = XX3, I_v = XX5)
DF5 = DataFrame(t = time,S_p = YY1, I_p = YY3, I_v = YY5)

CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Persistence_Q50p_low.csv",DF4)
CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1//Persistence_Q50p_high.csv",DF5)
=#
