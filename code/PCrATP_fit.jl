using Markdown, JSON
using DifferentialEquations, ModelingToolkit, Sundials, Interpolations
using DataFrames, VegaLite
using Revise
includet("../oxphos/OxPhos.jl")
using .OxPhos
include("../../simulation_functions.jl")
include("../../vl_functions.jl")

##
md"""
Load in data for PCr/ATP ratio
"""
(VO2_PCrATP,PCrATP_data) = OxPhos.load_PCrATP_data()

##
md"""
Plot the data
"""
fig = @vlplot(config={
        view={stroke=:grey},
        legend={title=""},
        font="Arial"
    },
    {:point, filled=true}, width=200, height=150,
    x={
        VO2_PCrATP, 
        axis={title="VO2 (μmol/min/g dry wt)", titleFontWeight=:normal, grid=false, tickSize=3}
    }, 
    y={
        PCrATP_data, 
        axis={title="PCr/ATP", titleFontWeight=:normal, grid=false, tickSize=3}
    }
)

enlarge(fig)

##
md"""
Fit the model's μa_ANT parameter to PCr/ATP ratio at different workloads
"""
popt,(μa_ANT_vals,costs) = OxPhos.fit_ANT_to_PCrATP()

##
md"""
Plot cost against parameter
"""
sys = bg_model_invivo()
copt = OxPhos.PCr_cost(sys,popt)
df = DataFrame(p=μa_ANT_vals,c=costs)
fig = df |> @vlplot(
    config={
        view={stroke=:grey},
        legend={title=""},
        font="Arial"
    }
) + 
@vlplot(
    :line, width=200, height=150,
    x={
        :p,
        axis={title="μa_ANT (kJ/mol)", titleFontWeight=:normal, grid=false, tickSize=3}
    }, 
    y={
        :c,
        axis={title="Cost", titleFontWeight=:normal, grid=false, tickSize=3}
    }
) + 
@vlplot({:point, color=:black, filled=true},x={datum=popt},y={datum=copt})
enlarge(fig)

##
md"""
Plot results for optimised parameter
"""
VO2, PCrATP = OxPhos.simulate_VO2(sys,popt)
df1 = DataFrame(VO2=VO2_PCrATP,PCrATP=PCrATP_data,type=:data)
df2 = DataFrame(VO2=VO2,PCrATP=PCrATP,type=:model)
df = [df1;df2]

fig = df |> @vlplot(
    config={
        view={stroke=:grey},
        legend={title=""},
        font="Arial"
    }
) + 
@vlplot(
    {:point, filled=true}, width=200, height=150,
    transform=[
        {filter="datum.type == 'data'"}
    ],
    x={
        :VO2, 
        axis={title="VO2 (μmol/min/g dry wt)", titleFontWeight=:normal, grid=false, tickSize=3}
    }, 
    y={
        :PCrATP, 
        axis={title="PCr/ATP ratio", titleFontWeight=:normal, grid=false, tickSize=3}
    }
) + 
@vlplot(
    mark={:line, color=:black},
    transform=[
        {filter="datum.type == 'model'"},
        {filter="datum.VO2 < 150"}
    ],
    x=:VO2, 
    y=:PCrATP
)

savevl(fig,"figures/PCrATP_fit")
enlarge(fig)
