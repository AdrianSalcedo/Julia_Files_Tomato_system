using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using DataFrames
using StochasticDiffEq

beta_p = 0.1
r_1 = 0.7
r_2 = 0.4
b = 0.5
beta_v = 0.27718623
theta = 0.818469652
mu = 0.058596601
gamma = 0.6
sigma_L = 0.378173863
sigma_I = 0.030015876
sigma_v = 0.254146903
N_v = mu/gamma
r= max(r_1,r_2)
u_0 = [70.0,20.0,10.0,3.0,4.0]
T = 200.0
t_s = range(0.0,T, step=1.0)
time = (0.0,T)
N_p = u_0[1]+u_0[2]+u_0[3]
dt=0.001

Rs0 = beta_p*beta_v/(r*gamma)
print("Rs0=",Rs0)

function F_Det(du,u,p,t)
 @inbounds begin
        du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
        du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
        du[3] = b*u[2]-r_2*u[3]
        du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
        du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
    end
    nothing
end

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
################################################################################
######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Det,u_0,time)
det_sol = solve(prob_det,Tsit5())
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,noise_rate_prototype=zeros(5,2))
sol = solve(prob_sde_tomato_sys,EM(),dt= dt)
################################################################################
############################ PLot variables ####################################
################################################################################

title = plot(title = "R_s =$Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
p1=plot(det_sol,vars=(1),color="blue")
p1=plot!(sol,vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(det_sol,vars=(2),color="blue")
p2=plot!(sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(det_sol,vars=(3),color="blue")
p3=plot!(sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(det_sol,vars=(4),color="blue")
p4=plot!(sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(det_sol,vars=(5),color="blue")
p5=plot!(sol,vars=(5),color="red",title ="Infec. v.")

plot(p1,p2,p3,p4,p5,title,layout = @layout([ [A B C]; [D E F]]), label="")

################################################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################
Datos=DataFrame()

j = 0
trajectories = 1
while j <= 10000
    monte_prob = MonteCarloProblem(prob_sde_tomato_sys)
    sim = solve(monte_prob, SROCKC2(),dt= dt,EnsembleThreads(),trajectories=trajectories)
    component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
    component = transpose(component) #transpose to obtain any*5 data matrix
    component = vcat(component...) #to obtain shape for dataframe
    component = vcat(component...) # again do a reshape
    variables = DataFrame(component) # define first data frame
    Datos_aux = DataFrame(t = t_s, S_p = variables[:,1], I_p = variables[:,3], I_v = variables[:,5]) #only some variables
    Datos = append!(Datos, Datos_aux)#append the data in the loop
    j+=1
    println("acepted =",j)
end

#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/R0+1/Trajectories//Data11.csv",Datos2)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/R0-1/Trajectories//Data.csv",Datos)
