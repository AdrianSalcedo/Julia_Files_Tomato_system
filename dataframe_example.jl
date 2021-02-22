using DataFrames
using CSV
df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
#DataFrame(CSV.File(input))
#df = DataFrame(x = 1, y = 2)
CSV.write("C:/Users/adria/Dropbox/Artículos/JuliaPro_code//Example.csv",df)
