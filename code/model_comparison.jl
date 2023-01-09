using ModelingToolkit, DifferentialEquations, VegaLite, DataFrames, Sundials
include("OxPhos.jl")
using .OxPhos
include("simulation_functions.jl")
include("vl_functions.jl")

## Compare a bond graph model with original model
@parameters t
alg = CVODE_BDF()

sys_beard = beard_model()
prob_beard = ODEProblem(sys_beard,[],(0.0,60.0))
sol_beard = solve(prob_beard,alg; reltol=1e-7, abstol=1e-10)

sys_bg = bg_model_invitro()
prob_bg = ODEProblem(sys_bg,[],(0.0,60.0))
sol_bg = solve(prob_bg,alg; reltol=1e-7, abstol=1e-10)

@variables dPsi(t)
df1 = DataFrame(t=sol_beard.t,dPsi=sol_beard[dPsi],model="Kinetic")
df2 = DataFrame(t=sol_bg.t,dPsi=sol_bg[dPsi],model="Bond graph")
df = [df1;df2]
fig = df |> @vlplot(
    :line, width=240, height=180,
    x={
        :t, 
        axis={title="Time (s)", titleFontWeight=:normal, grid=false}}, 
    y={
        :dPsi, 
        axis={
            title="dPsi (mV)", titleFontWeight=:normal, grid=false, offset=5},
        scale={zero=false}
    },
    color={:model, legend={title=""}, scale={scheme=:set2}},
    config={view={stroke=:transparent}}
)
enlarge(fig)

## Plot comparison at different external Pi concentrations
@parameters ADP_e Pi_e
Pi_concentrations=1e-3*(0:0.1:10)
ps = [[ADP_e => 1.3e-3, Pi_e => p] for p in Pi_concentrations]

tspan = (0.0,60.0)
sols = multi_sim(sys_beard,tspan,ps)
dPsi_beard = [sol[dPsi][end] for sol in sols]

sols = multi_sim(sys_bg,tspan,ps)
dPsi_bg = [sol[dPsi][end] for sol in sols]


df1 = DataFrame(Pi=1000*Pi_concentrations,dPsi=dPsi_beard,model="Kinetic")
df2 = DataFrame(Pi=1000*Pi_concentrations,dPsi=dPsi_bg,model="Bond graph")
df = [df1;df2]
fig = df |> @vlplot(
    {:line, color=:blue}, width=240, height=180,
    x={
        :Pi, 
        axis={title="[Pi]ₑ (mM)", titleFontWeight=:normal, grid =false}}, 
    y={
        :dPsi, 
        axis={
            title="dPsi (mV)", titleFontWeight=:normal, grid =false},
        scale={zero=false}
    },
    color={:model, legend={title=""}, scale={scheme=:set2}},
    config={view={stroke=:transparent}}
)
enlarge(fig)



