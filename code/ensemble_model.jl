using ModelingToolkit, DifferentialEquations, Sundials, Markdown, Random, JLD2, Setfield, Interpolations
include("simulation_functions.jl")
include("vl_functions.jl")
include("OxPhos.jl")
include("data_distributions.jl")
using .DataDistributions

Random.seed!(2702)
run_simulations = false
output_path = "../output/ensemble_simulations"

## 
md"""
Load in fold changes
"""
include("ageing_data.jl")
omics_path = "../data/omics4path.csv"
oxphos_protein_path = "../gene_names/oxphos_proteins.xlsx"
oxphos_metab_path = "../gene_names/oxphos_metabolites.xlsx"
(df_reaction_FC, df_metabolite_FC) = extract_FCs(omics_path, oxphos_protein_path, oxphos_metab_path)
proteomics_info = extract_proteomics(omics_path,oxphos_protein_path)

## 
md"""
Extract variances of metabolites
"""
metab_data = read_data("../data/metabAll_log2EigenMSImp.csv")
metab_names = ["Nicotinamide adenine dinucleotide", "Reduced nicotinamide Adenine Dinucleotide"]
metab_data_filter = filter_names(metab_data,metab_names)
(young_metab_data,old_metab_data) = DataDistributions.split(metab_data_filter,12)


## 
md"""
Extract variances of proteins
"""
prot_data = read_data("../data/proteins_log2Imp.csv")
prot_names = ["SLC25A4", "LETM1", "SLC25A3", "CKM", "CKMT2"]
prot_data_filter = filter_names(prot_data,prot_names)
(young_prot_data,old_prot_data) = DataDistributions.split(prot_data_filter,12)

##
md"""
Apply the simulated variations to the parameter values and run simulations
"""
n_sims = 10000
young_metab_ratios = sim2Gaussian(young_metab_data, n_sims, power=true)
old_metab_ratios = sim2Gaussian(old_metab_data, n_sims, power=true)
young_prot_ratios = sim2Gaussian(young_prot_data, n_sims, power=false)
old_prot_ratios = sim2Gaussian(old_prot_data, n_sims, power=false)

@parameters t
@variables NAD_x(t) NADH_x(t)

include("simulations.jl")
sys_y = ageing_sys()
sys_o = ageing_sys()
u0_o = adjust_u0()
sys_o.defaults[NAD_x] = u0_o[1][2]
sys_o.defaults[NADH_x] = u0_o[2][2]

# Run ensemble simulations
metabolite_map = Dict(
    "Nicotinamide adenine dinucleotide" => NAD_x,
    "Reduced nicotinamide Adenine Dinucleotide" => NADH_x
)
@parameters μa_ANT μa_KH μa_PiHt μa_CK μa_CKi
protein_map = Dict(
    "SLC25A4" => μa_ANT, 
    "LETM1" => μa_KH, 
    "SLC25A3" => μa_PiHt, 
    "CKM" => μa_CK, 
    "CKMT2" => μa_CKi
)

function generate_inputs(sys,index,metab_matrix,prot_matrix)
    u0 = Pair{Num, Float64}[]
    for (i,metab) in enumerate(metab_names)
        variable = metabolite_map[metab]
        ratio = metab_matrix[i,index]
        push!(u0, variable => ratio*sys.defaults[variable])
    end

    p = Pair{Num, Float64}[]
    for (i,prot) in enumerate(prot_names)
        parameter = protein_map[prot]
        # Since the variables are drawn from log2-space and we want thermodynamic parameters in ln-space, we need to convert by a factor of RT*ln(2)
        change = prot_matrix[i,index] * OxPhos.RT_default * log(2)
        push!(p, parameter => sys.defaults[parameter] + change)
    end
    ps = [[p_ATPase; p] for p_ATPase in workload_param_vec(n_steps=50,max_wl=20)]

    return (u0,p,ps)
end

function simulate_variation(sys,index,metab_matrix,prot_matrix)
    (u0,p,ps) = generate_inputs(sys,index,metab_matrix,prot_matrix)

    tspan = (0.0,1000.0)
    options = Dict(:reltol => 1e-10, :abstol => 1e-15, :alg => CVODE_BDF())
    return multi_sim(sys,tspan,ps; u0=u0,options...)
end

function extract_vars(sols,group,index)
    @variables PCr(t) Cr(t) ATP_fe_d(t) ATP_me_d(t) J_C4(t)
    VO2 = 0.5*extract_var(sols,J_C4)
    Cr_e = extract_var(sols,Cr)
    PCr_e = extract_var(sols,PCr)

    @variables ADP_fe_d(t) ADP_me_d(t) ATP_fe_d(t) ATP_me_d(t) J_C4(t) μ_ATP_me(t) μ_ADP_me(t) μ_Pi_e(t) μ_H_i(t)
    @parameters μ_H2O
    ADP = extract_var(sols,ADP_fe_d+ADP_me_d)
    ATP = extract_var(sols,ATP_fe_d+ATP_me_d)
    
    return DataFrame(
        group=group,
        index=index,
        VO2=VO2, 
        VO2_m=1000*VO2, 
        PCrATP=PCr_e./ATP,
        Cr=Cr_e,
        PCr=PCr_e,
        ATP=ATP,
        ADP=ADP,
        ADP_μ=1e6*ADP,
        ADPATP=ADP./ATP,
    )
end

function run_and_save(sys,metab_ratios,prot_ratios,path,group)
    for index in 1:n_sims
        sols = simulate_variation(sys,index,metab_ratios,prot_ratios)
        df = extract_vars(sols,group,index)
        jldopen("$path/$(string(group))_sols.jld2", "a+") do file
            file[string(index)] = [sol[end] for sol in sols]
        end
        jldopen("$path/$(string(group))_dfs.jld2", "a+") do file
            file[string(index)] = df
        end
    end
end

if run_simulations
    files = ["Young_sols","Young_dfs","Old_sols","Old_dfs"]
    for file in files
        jldopen("$output_path/$file.jld2", "w") do file
        end
    end

    run_and_save(sys_y,young_metab_ratios,young_prot_ratios,output_path,:Young)
    run_and_save(sys_o,old_metab_ratios,old_prot_ratios,output_path,:Old)
end

dfs_y = [load("$output_path/Young_dfs.jld2",string(i)) for i in 1:n_sims]
dfs_o = [load("$output_path/Old_dfs.jld2",string(i)) for i in 1:n_sims]
df = vcat(dfs_y...,dfs_o...)
    
##
md"""
Plot the results
"""
line_plot(x,y) = @vlplot(
    config={
        view={stroke=:transparent},
        legend={title="",orient="top-right",offset=5},
        font="Arial"
    }
) + @vlplot(
    mark={:line, strokeWidth=0.5, opacity=0.5}, 
    width=120, height=120,
    x={
        x, 
        axis={
            titleFontWeight=:normal,
            grid=false,
            tickSize=3,
            domainColor=:black,
            tickColor=:black
        },
        scale={zero=true}
    }, 
    y={
        y, 
        axis={
            titleFontWeight=:normal, 
            grid=false,
            tickSize=3,
            domainColor=:black,
            tickColor=:black
        },
        scale={zero=true}
    },
    detail = :index,
    color={:group, scale={range=["#009DA3","#F98F00"]}, sort=:descending}
)
include("figure_adjustments.jl")

df_plot = subset(df, :index => ByRow(i -> i <= 200))

fig = df_plot |> line_plot(:VO2_m,:PCrATP)
set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig,"y","PCr:ATP ratio")
savevl(fig,"$output_path/ensemble_PCrATP")

fig = df_plot |> line_plot(:VO2_m,:ADP_μ)
set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig,"y","ADP (μM)")
set_domain!(fig,"x",[0,3])
set_domain!(fig,"y",[0,100])
fig = @set fig.config.legend.orient = "bottom-right"
savevl(fig,"$output_path/ensemble_ADP")

##
"""
Average curves over the samples generated for each group
"""
function check(dfs,i)
    VO2 = dfs[i].VO2
    return all(isnan.(VO2) .== false)
end
extract_mean(dfs,var,j) = mean(dfs[i][!,var][j] for i in eachindex(dfs))
extract_means(dfs,var) = [extract_mean(dfs,var,j) for j in eachindex(dfs[1].VO2)]

function summarise_means(dfs,group)
    mean_VO2 = extract_means(dfs,:VO2_m)
    mean_PCrATP = extract_means(dfs,:PCrATP)
    mean_ADP = extract_means(dfs,:ADP_μ)
    return DataFrame(VO2=mean_VO2,PCrATP=mean_PCrATP,ADP=mean_ADP,group=group)
end

df_mean_y = summarise_means(dfs_y,:Young)
df_mean_o = summarise_means(dfs_o,:Old)
df_mean = vcat(df_mean_y,df_mean_o)

fig_PCrATP = df_mean |> line_plot(:VO2,:PCrATP)
fig_PCrATP.layer[1]["mark"]["opacity"] = 1
delete!(fig_PCrATP.layer[1]["mark"],"strokeWidth")
set_axis_title!(fig_PCrATP,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig_PCrATP,"y","PCr:ATP ratio")
savevl(fig_PCrATP,"$output_path/ensemble_mean_PCrATP")

fig_ADP = df_mean |> line_plot(:VO2,:ADP)
fig_ADP.layer[1]["mark"]["opacity"] = 1
delete!(fig_ADP.layer[1]["mark"],"strokeWidth")
set_axis_title!(fig_ADP,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig_ADP,"y","ADP (μM)")
set_domain!(fig_ADP,"y",[0,300])
fig_ADP = @set fig_ADP.config.legend.orient = "bottom-right"
savevl(fig_ADP,"$output_path/ensemble_mean_ADP")

##
"""
Plot distributions of metabolites
"""
# Load in final simulations states
u_end_y = [load("$output_path/Young_sols.jld2",string(i)) for i in 1:n_sims]
u_end_o = [load("$output_path/Old_sols.jld2",string(i)) for i in 1:n_sims]

# Try another approach... I should have all the data I need already
all_metabolites = @variables NAD_x(t) NADH_x(t) Q(t) QH2(t) Cred(t) Cox(t) Pi_x(t) Pi_i(t) Pi_e_d(t) H_x(t) Mg_x(t) Mg_e_d(t) K_x(t) Cr(t) PCr(t) Cr_i(t) PCr_i(t) ATP_fx(t) ATP_mx(t) ATP_fi(t) ATP_mi(t) ATP_fe_d(t) ATP_me_d(t) ADP_fx(t) ADP_mx(t) ADP_fi(t) ADP_mi(t) ADP_fe_d(t) ADP_me_d(t)
indexof(sym,syms) = findfirst(isequal(sym),syms)
metabolite_index_map = Dict(
    x => indexof(x,states(sys_y)) for x in all_metabolites
)

# Define a function to extract the right VO2
VO2_human = 6.05e-4

splice(x,index) = [x[i][index] for i in 1:length(x)]

@variables Q(t) QH2(t) Cox(t) Cred(t)
@parameters Qtot Ctot
pb = OxPhos.parameters_beard()
qtot = pb[Qtot]
ctot = pb[Ctot]

function interpolate(x_vec,VO2_vec,VO2)
    interp = LinearInterpolation(VO2_vec,x_vec)
    x = interp(VO2)
end

process_name(m::Num) = string(m)[1:end-3]
process_name(m::String) = m

function insert_metabolite!(dict,m,x,group)
    push!(dict[:metabolite], process_name(m))
    push!(dict[:concentration], x)
    push!(dict[:group], group)
    return
end

function interpolate_and_insert!(dict,m,group,x_vec,VO2_vec,VO2)
    x = 0
    try
        x = interpolate(x_vec,VO2_vec,VO2)
    catch
        x = NaN
    end
    insert_metabolite!(dict,m,x,group)
end

function sum_metabolite_concentrations(array_u,ms,metabolite_index_map)
    return sum(splice(array_u,metabolite_index_map[m]) for m in ms)
end

function update_metabolite_concentrations!(dict,df,array_u,group;
    metabolite_index_map=metabolite_index_map,VO2=2*VO2_human)
    VO2_vec = df.VO2

    # Handle single metabolites
    x = 0

    single_metabolites = (
        @variables NAD_x(t) NADH_x(t) Q(t) QH2(t) Cred(t) Cox(t) Pi_x(t) Pi_i(t) Pi_e_d(t) H_x(t) Mg_x(t) Mg_e_d(t) K_x(t) Cr(t) PCr(t) Cr_i(t) PCr_i(t)
    )

    for m in single_metabolites
        i = metabolite_index_map[m]
        if i isa Number
            x_vec = splice(array_u,i)
        elseif m === Q
            x_vec = qtot .- splice(array_u,metabolite_index_map[QH2])
        elseif m === Cox
            x_vec = ctot .- splice(array_u,metabolite_index_map[Cred])
        end
        interpolate_and_insert!(dict,m,group,x_vec,VO2_vec,VO2)
    end

    # Handle lumped metabolites
    @variables ATP_fx(t) ATP_mx(t) ATP_fi(t) ATP_mi(t) ATP_fe_d(t) ATP_me_d(t) ADP_fx(t) ADP_mx(t) ADP_fi(t) ADP_mi(t) ADP_fe_d(t) ADP_me_d(t)
    # ATP_x
    x_vec = sum_metabolite_concentrations(array_u,(ATP_fx, ATP_mx),metabolite_index_map)
    interpolate_and_insert!(dict,"ATP_x",group,x_vec,VO2_vec,VO2)
    
    # ATP_i
    x_vec = sum_metabolite_concentrations(array_u,(ATP_fi, ATP_mi),metabolite_index_map)
    interpolate_and_insert!(dict,"ATP_i",group,x_vec,VO2_vec,VO2)
    #x = interpolate(x_vec,VO2_vec,VO2)
    #insert_metabolite!(dict,"ATP_i",x,group)

    # ATP_e
    x_vec = sum_metabolite_concentrations(array_u,(ATP_fe_d, ATP_me_d),metabolite_index_map)
    interpolate_and_insert!(dict,"ATP_e",group,x_vec,VO2_vec,VO2)
    #x = interpolate(x_vec,VO2_vec,VO2)
    #insert_metabolite!(dict,"ATP_e",x,group)

    # ADP_x
    x_vec = sum_metabolite_concentrations(array_u,(ADP_fx, ADP_mx),metabolite_index_map)
    interpolate_and_insert!(dict,"ADP_x",group,x_vec,VO2_vec,VO2)
    #x = interpolate(x_vec,VO2_vec,VO2)
    #insert_metabolite!(dict,"ADP_x",x,group)

    # ADP_i
    x_vec = sum_metabolite_concentrations(array_u,(ADP_fi, ADP_mi),metabolite_index_map)
    interpolate_and_insert!(dict,"ADP_i",group,x_vec,VO2_vec,VO2)
    #x = interpolate(x_vec,VO2_vec,VO2)
    #insert_metabolite!(dict,"ADP_i",x,group)

    # ADP_e
    x_vec = sum_metabolite_concentrations(array_u,(ADP_fe_d, ADP_me_d),metabolite_index_map)
    interpolate_and_insert!(dict,"ADP_e",group,x_vec,VO2_vec,VO2)
    #x = interpolate(x_vec,VO2_vec,VO2)
    #insert_metabolite!(dict,"ADP_e",x,group)

    return dict
end

dict_metabolites = OrderedDict(:metabolite => String[], :group => Symbol[], :concentration => Float64[])

for (df,array_u) in zip(dfs_y,u_end_y)
    update_metabolite_concentrations!(dict_metabolites,df,array_u,:Young)
end
for (df,array_u) in zip(dfs_o,u_end_o)
    update_metabolite_concentrations!(dict_metabolites,df,array_u,:Old)
end

df_metabolites = DataFrame(dict_metabolites)
df_metabolites.concentration .*= 1000

function plot_distribution(df_metabolites,metabolite)
    df_subset = subset(subset(df_metabolites, :metabolite => ByRow(i -> i == metabolite)), :concentration => ByRow(x -> isnan(x) == false))

    x_min = quantile(df_subset.concentration, 0.02)
    x_max = quantile(df_subset.concentration, 0.98)
    bandwidth = maximum(x_max-x_min)/100

    df_y = subset(df_subset, :group => ByRow(g -> g == :Young))
    df_o = subset(df_subset, :group => ByRow(g -> g == :Old))
    median_y = median(df_y.concentration)
    median_o = median(df_o.concentration)

    fig = df_subset |> @vlplot(
        config={
            view={stroke=:transparent},
            legend={title="",orient="top-right",offset=5},
            font="Arial"
        }
    )+ @vlplot(
        mark={:area, opacity=0.5, clip=true},
        width=160,
        height=120,
        transform=[
            {
                density=:concentration,
                bandwidth=bandwidth,
                extent=[x_min,x_max],
                groupby=[:group]
            }
        ],
        x={
            "value:q",
            axis={
                title="$metabolite (mM)",
                titleFontWeight=:normal, 
                grid=false,
                tickSize=3,
                domainColor=:black,
                tickColor=:black
            },
            scale={
                domain=[x_min,x_max]
            }
        },
        y={
            "density:q",
            axis={
                title="Density",
                titleFontWeight=:normal, 
                grid=false,
                tickSize=3,
                domainColor=:black,
                tickColor=:black
            }
        },
        color={"group:n", scale={range=["#009DA3","#F98F00"]}, sort=:descending}
    ) + 
    @vlplot({:rule, color="#009DA3"}, x={datum=median_y}) +
    @vlplot({:rule, color="#F98F00"}, x={datum=median_o})

    if metabolite == "H_x"
        fig.params["layer"][1]["encoding"]["x"]["axis"]["format"] = ".0e"
    end

    return fig
end

function plot_box(df_metabolites,metabolite)
    df_subset = subset(subset(df_metabolites, :metabolite => ByRow(i -> i == metabolite)), :concentration => ByRow(x -> isnan(x) == false))

    fig = df_subset |> @vlplot(
        config={
            view={stroke=:transparent},
            legend={title="",orient="top-right",offset=5},
            font="Arial"
        }
    )+ @vlplot(
        mark={:boxplot, extend=1.5, orient=:horizontal},
        title=metabolite,
        width=120,
        height=60,
        y={
            "group:n",
            axis=nothing,
            sort=:descending
        },
        x={
            "concentration:q",
            axis={
                titleFontWeight=:normal, 
                grid=false,
                tickSize=3,
                domainColor=:black,
                tickColor=:black
            }
        },
        color={"group:n", scale={range=["#009DA3","#F98F00"]}, sort=:descending}
    )
end

function plot_metabolomics_results(df_metabolites)
    figs = Dict(
        m => plot_distribution(df_metabolites,m) for m in unique(df_metabolites.metabolite)
    )
    
    collated_fig = @vlplot(config={
        view={stroke=:transparent},
        legend={title=""},
        font="Arial"
    }) + [
        [figs["NAD_x"] figs["NADH_x"] figs["Q"] figs["QH2"]];
        [figs["Cred"] figs["Cox"] figs["Pi_x"] figs["Pi_i"]];
        [figs["Pi_e_d"] figs["H_x"] figs["Mg_x"] figs["Mg_e_d"]];
        [figs["K_x"] figs["Cr"] figs["PCr"] figs["Cr_i"]];
        [figs["PCr_i"] figs["ADP_e"] figs["ATP_x"] figs["ATP_e"]];
        [figs["ADP_i"] figs["ATP_i"] figs["ADP_x"]]
    ]
    return collated_fig
end

fig_ensemble = plot_metabolomics_results(df_metabolites)
savevl(fig_ensemble,"$output_path/metabolites_high_wl")