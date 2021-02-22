using Plots; plotly()
using CSV
using IterableTables, DataTables
using DataFrames
using Statistics

Data = DataFrame()
#Windows version
#=for i in 1:100
    #Data_aux = CSV.read("\\home\\gabrielsalcedo\\Dropbox\\Artículos\\JuliaPro_code\\Noise_Extinction\\Trajectories\\Data$i.csv")
    Data_aux = CSV.read("C:\\Users\\adria\\Dropbox\\Artículos\\JuliaPro_code\\Noise_Extinction\\Trajectories\\Data$i.csv")
    Data = append!(Data, Data_aux)
end
=#
#Linux version
#Data = CSV.read("C:\\Users\\adria\\Dropbox\\Artículos\\JuliaPro_code\\Noise_Extinction\\Trajectories\\Data_new_2.csv",DataFrame)
#Data = CSV.read("\\home\\gabrielsalcedo\\Dropbox\\Artículos\\JuliaPro_code\\Noise_Extinction\\Trajectories\\Data02.csv", DataFrame)
#CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/Noise_Extinction/Trajectories//Data.csv",Data)
Data = CSV.read("D:\\Data_new_4.csv",DataFrame)


Data_Mean = DataFrame()

for i in 0:200
    Sub = Data[Data.t .== i, :]
    Mean_t = mean(Sub.t)
    Mean_S_p = mean(Sub.S_p)
    Mean_I_p = mean(Sub.I_p)
    Mean_I_v = mean(Sub.I_v)
    Quartile_S_p = quantile!(Sub.S_p, [0.05, 0.25, 0.75, 0.95])
    Quartile_I_p = quantile!(Sub.I_p, [0.05, 0.25, 0.75, 0.95])
    Quartile_I_v = quantile!(Sub.I_v, [0.05, 0.25, 0.75, 0.95])
    Data_aux = DataFrame(
    t = Mean_t,
    Q05_S_p = Quartile_S_p[1],
    Q25_S_p = Quartile_S_p[2],
    Mean_S_p = Mean_S_p,
    Q75_S_p = Quartile_S_p[3],
    Q95_S_p = Quartile_S_p[4],
    Q05_I_p = Quartile_I_p[1],
    Q25_I_p = Quartile_I_p[2],
    Mean_I_p = Mean_I_p,
    Q75_I_p = Quartile_I_p[3],
    Q95_I_p = Quartile_I_p[4],
    Q05_I_v = Quartile_I_v[1],
    Q25_I_v = Quartile_I_v[2],
    Mean_I_v = Mean_I_v,
    Q75_I_v = Quartile_I_v[3],
    Q95_I_v = Quartile_I_v[4],
      )
    Data_Mean = append!(Data_Mean,Data_aux)
end
CSV.write("D://Data_Mean_Noise_Exctinction.csv",Data_Mean)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code/Noise_Extinction/Trajectories//Data_Mean_Extinction_1.csv",Data_Mean)
#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/Noise_Extinction/Trajectories//Data_Mean_Extinction.csv",Data_Mean)
#\\home\\gabrielsalcedo\\Dropbox\\Artículos\\JuliaPro_code\\Noise_Extinction\\Trajectories
