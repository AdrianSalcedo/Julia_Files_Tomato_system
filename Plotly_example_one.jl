using PlotlyJS
#Plotly.set_credentials_file(Dict("username"=>"AdrianSalcedo","api_key"=>"••••••••••"))

data = scatter(;
    x = ["2013-10-04 22:23:00", "2013-11-04 22:23:00", "2013-12-04 22:23:00"],
    y = [1, 3, 6],
    mode = "lines")

layout = Layout(;title="Data Labels Hover",
                 xaxis_range=[0, 10],
                 yaxis_range=[0, 7])

PlotlyJS.plot(data, layout)

response = PlotlyJS.plot(data, ["filename" => "date-axes", "fileopt" => "overwrite"])
plot_url = response["url"]
