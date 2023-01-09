using ModelingToolkit, DifferentialEquations, Sundials, Markdown, DataFrames, VegaLite
using CSV, XLSX, OrderedCollections, Cleaner, Statistics, Interpolations
using Plots
include("OxPhos.jl")
using .OxPhos
include("simulation_functions.jl")
include("vl_functions.jl")

## Model definition
@parameters t
@variables NAD_x(t) NADH_x(t) J_ATPase(t)
@parameters μa_ANT μa_ATPase

# Import the model
tspan = (0.0,1000.0)
sys = bg_model_invivo(stat=[NAD_x,NADH_x],species=:rat)
sys.defaults[μa_ATPase] = L(5)+sys.defaults[μa_ATPase]
options = Dict(:reltol => 1e-7, :abstol => 1e-10, :alg => CVODE_BDF())

## Do a sensitivity analysis without fitting (on the rat model), at a high workload
function sensitivity_sims(sys,tspan,p,ratios;u0=[])
    defaults = ModelingToolkit.get_defaults(sys)
    ps = [[p => L(r)+defaults[p]] for r in ratios]
    return multi_sim(sys,tspan,ps;options...)
end

function perturb_parameters(sys,ratios)
    perturbation_parameters = @parameters μa_C1 μa_C3 μa_C4 μa_F1 μa_MgATP μa_MgADP μa_Hle μa_KH μa_DH μa_fATPt μa_fADPt μa_mATPt μa_mADPt μa_AMPt μa_Pit μa_AK μa_ANT μa_PiHt μa_CK μa_CKi μ_r μa_ATPase
    sols = Dict(
        p => sensitivity_sims(sys,tspan,p,ratios) 
        for p in perturbation_parameters
    );
    df = DataFrame(parameter=Num[], ratio_down=Float64[], ratio_up=Float64[])
    for p in perturbation_parameters
        v = extract_var(sols[p],J_ATPase)
        v_base = v[length(v)÷2 + 1]
        push!(df,(p, v[1]/v_base-1, v[end]/v_base-1))
    end
    return sort(df, order(:ratio_down, by=abs, rev=true))
end

ratios = [0.7, 1.0, 1.3]
perturb_parameters(sys,ratios)

## Do a sensitivity analysis after fitting
# Find the ATPase parameter to use
@variables J_C4(t)
sys = bg_model_invivo(stat=[NAD_x,NADH_x])
sys.defaults[μa_ANT] = OxPhos.updated_μa_ANT()
VO2_human = 6.05e-4
VO2_sim = 2*VO2_human
# Interpolate to find the rates of ATPases
ATP_consumption_ratio = range(0.1,stop=7,length=100)
ps = OxPhos.parameters_bg_invivo()
ATPase_param_vals = @. ps[μa_ATPase] + L(ATP_consumption_ratio)
p = [[μa_ATPase => m] for m in ATPase_param_vals]
sols = multi_sim(sys,tspan,p; options...)
VO2_vec = [0.5*sol[J_C4][end] for sol in sols]
function find_ATPase_rate(VO2,VO2_vec,ATPase_rates)
    interp = LinearInterpolation(VO2_vec,ATPase_rates)
    return interp(VO2)
end
μa_ATPase_val = find_ATPase_rate(VO2_sim,VO2_vec,ATPase_param_vals)
sys.defaults[μa_ATPase] = μa_ATPase_val

ratios = [0.7, 1.0, 1.3]
perturb_parameters(sys,ratios)