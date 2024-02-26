using ModelingToolkit, DifferentialEquations, Sundials, Markdown, DataFrames, VegaLite
using Interpolations, Setfield
include("OxPhos.jl")
include("simulation_functions.jl")
include("vl_functions.jl")
include("figure_adjustments.jl")
output_vo2 = "../output/VO2"

## 
md"""
Load in data
"""
include("ageing_data.jl")
omics_path = "../data/omics/omics4path.csv"
oxphos_protein_path = "../gene_names/oxphos_proteins.xlsx"
oxphos_metab_path = "../gene_names/oxphos_metabolites.xlsx"
(df_reaction_FC, df_metabolite_FC) = extract_FCs(omics_path, oxphos_protein_path, oxphos_metab_path)

##
md"""
Run simulations of the model for young and old hearts
"""
include("simulations.jl")

array_FC_scale = [0.0,1.0]
sys = ageing_sys()
array_sols = [:Young => simulate_FCs(sys,0,0), :Old => simulate_FCs(sys,1,1)];

## 
"""
Extract variables from solutions
"""
function extract_vars(sols,group)
    @variables PCr(t) Cr(t) ATP_fe_d(t) ATP_me_d(t) J_C4(t)
    VO2 = 0.5*extract_var(sols,J_C4)
    Cr_e = extract_var(sols,Cr)
    PCr_e = extract_var(sols,PCr)

    @variables dPsi(t) μ_H_x(t) μ_H_i(t) H_x(t)
    F = 0.096484
    dPsi = extract_var(sols,dPsi)
    H_x = extract_var(sols,H_x)
    μ_H_x = extract_var(sols,μ_H_x)
    μ_H_i = extract_var(sols,μ_H_i)
    Δp = dPsi + (μ_H_i - μ_H_x)/F

    @variables ADP_fe_d(t) ADP_me_d(t) ATP_fe_d(t) ATP_me_d(t) J_C4(t) μ_ATP_me(t) μ_ADP_me(t) μ_Pi_e(t) μ_H_i(t)
    @parameters μ_H2O
    ADP = extract_var(sols,ADP_fe_d+ADP_me_d)
    ATP = extract_var(sols,ATP_fe_d+ATP_me_d)

    μ_ATP_me = extract_var(sols,μ_ATP_me)
    μ_H2O = OxPhos.bg_species_parameters()[μ_H2O]
    μ_ADP_me = extract_var(sols,μ_ADP_me)
    μ_Pi_e = extract_var(sols,μ_Pi_e)
    μ_H_i = extract_var(sols,μ_H_i)
    ΔG_ATP = @. μ_ADP_me + μ_Pi_e + μ_H_i - μ_ATP_me - μ_H2O
    
    return DataFrame(
        group=group,
        VO2=VO2, 
        VO2_m=1000*VO2, 
        PCrATP=PCr_e./ATP,
        Cr=Cr_e,
        PCr=PCr_e,
        ATP=ATP,
        ADP=ADP,
        ADP_μ=1e6*ADP,
        ADPATP=ADP./ATP,
        dPsi=dPsi,
        Hx=H_x,
        Δp=Δp,
        ΔG_ATP=ΔG_ATP
    )

end
dfs = [extract_vars(sols,string(scale)) for (scale,sols) in array_sols]
df = vcat(dfs...)

## 
"""
Plot PCr/ATP ratio against workload
"""
vl_config() = @vlplot(
    config={
        view={stroke=:transparent},
        legend={title="",orient="top-right",offset=5},
        font="Arial"
    }
)
vl_linespec(x,y) = @vlplot(
    :line, width=120, height=120,
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
    color={:group, scale={range=["#009DA3","#F98F00"]}, sort=:descending}
)
vl_line(df,x,y;config=vl_config()) = df |> config + vl_linespec(x,y)

function PCrATP_plot(df)
    fig = vl_line(df,:VO2_m,:PCrATP)
    set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
    set_axis_title!(fig,"y","PCr:ATP ratio")
    return fig
end

fig = PCrATP_plot(df)
savevl(fig,"$output_vo2/PCrATP")
enlarge(fig)

## 
"""
Plot membrane potential against workload
"""
fig = vl_line(df,:VO2_m,:dPsi)
set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig,"y","ΔΨ (mV)")
zero_axis!(fig,"y",false)

savevl(fig,"$output_vo2/dPsi")
enlarge(fig)

# Plot proton-motive force
fig = vl_line(df,:VO2_m,:Δp)
set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig,"y","Δp (mV)")
zero_axis!(fig,"y",false)

savevl(fig,"$output_vo2/dp")
enlarge(fig)

## 
"""
Plot ADP/ATP ratio (cytosolic) against workload
"""
fig = vl_line(df,:VO2_m,:ADPATP)
set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig,"y","ADP:ATP ratio")
set_domain!(fig,"y",[0,0.02])
fig = @set fig.config.legend.orient = "top-left"

savevl(fig,"$output_vo2/ADPATP")
enlarge(fig)

## 
"""
Plot ADP concentration (cytosolic) against workload
"""
function ADP_plot(df)
    fig = vl_line(df,:VO2_m,:ADP_μ)
    set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
    set_axis_title!(fig,"y","ADP (μM)")
    set_domain!(fig,"y",[0,100])
    fig = @set fig.config.legend.orient = "top-left"
    return fig
end

fig = ADP_plot(df)
savevl(fig,"$output_vo2/ADP")
enlarge(fig)

## 
"""
Plot free energy of ATP hydrolysis (cytosolic) against workload
"""
fig = vl_line(df,:VO2_m,:ΔG_ATP)
set_axis_title!(fig,"x","VO2 (mmol/L mito/s)")
set_axis_title!(fig,"y","ΔG_ATP (kJ/mol)")
zero_axis!(fig,"y",false)
fig = @set fig.config.legend.orient = "top-left"

savevl(fig,"$output_vo2/ATP_hydrolysis")
enlarge(fig)

##
md"""
Define functions for running simulations and plotting results
"""
# Interpolate to find the rates of ATPases
function find_ATPase_rate(VO2,VO2_vec,ATPase_rates)
    interp = LinearInterpolation(VO2_vec,ATPase_rates)
    return interp(VO2)
end

function simulate_metabolomics(VO2_sim)
    ATPase_param_vals = ATPase_params()
    μa_ATPase_y = find_ATPase_rate(VO2_sim,dfs[1].VO2,ATPase_param_vals) # Young
    μa_ATPase_o = find_ATPase_rate(VO2_sim,dfs[2].VO2,ATPase_param_vals) # Old

    prob_y = ODEProblem(sys,[],tspan,[μa_ATPase => μa_ATPase_y])
    sol_young = solve(prob_y;options...)

    prob_o = ODEProblem(sys,adjust_u0(1,1),tspan,[μa_ATPase => μa_ATPase_o])
    sol_old = solve(prob_o;options...)

    dict_fc = OrderedDict(:metabolite => String[], :young => Float64[], :old => Float64[], :FC => Float64[])
    single_metabolites = (
        @variables NAD_x(t) NADH_x(t) Q(t) QH2(t) Cred(t) Cox(t) Pi_x(t) Pi_i(t) Pi_e_d(t) H_x(t) Mg_x(t) Mg_e_d(t) K_x(t) Cr(t) PCr(t) Cr_i(t) PCr_i(t)
    )
    for m in single_metabolites
        x_young = sol_young[m][end]
        x_old = sol_old[m][end]
        FC = x_old/x_young

        push!(dict_fc[:metabolite], string(m)[1:end-3])
        push!(dict_fc[:young], x_young)
        push!(dict_fc[:old], x_old)
        push!(dict_fc[:FC], FC)
    end

    @variables ATP_fx(t) ATP_mx(t) ATP_fi(t) ATP_mi(t) ATP_fe_d(t) ATP_me_d(t) ADP_fx(t) ADP_mx(t) ADP_fi(t) ADP_mi(t) ADP_fe_d(t) ADP_me_d(t)
    lumped_metabolites = Dict(
        "ATP_x" => ATP_fx + ATP_mx,
        "ATP_i" => ATP_fi + ATP_mi,
        "ATP_e" => ATP_fe_d + ATP_me_d,
        "ADP_x" => ADP_fx + ADP_mx,
        "ADP_i" => ADP_fi + ADP_mi,
        "ADP_e" => ADP_fe_d + ADP_me_d
    )
    for (s,x) in lumped_metabolites
        x_young = sol_young[x][end]
        x_old = sol_old[x][end]
        FC = x_old/x_young

        push!(dict_fc[:metabolite], s)
        push!(dict_fc[:young], x_young)
        push!(dict_fc[:old], x_old)
        push!(dict_fc[:FC], FC)
    end

    df_metab_fc = DataFrame(dict_fc)
end

function plot_metabolomics_results(df_results)
    df_metab_fc_plot = DataFrame(df_results)
    df_metab_fc_plot.young .= 1000*df_metab_fc_plot.young
    df_metab_fc_plot.old .= 1000*df_metab_fc_plot.old
    df_metab_vars = stack(df_metab_fc_plot,[:young,:old],:metabolite)
    DataFrames.rename!(df_metab_vars,[:variable => :group, :value => :concentration])

    figs = Dict(
        m => bar_plot(m,df_metab_vars) for m in df_metab_vars.metabolite
    )
    
    collated_fig = @vlplot(config={
        view={stroke=:transparent},
        legend={title=""},
        font="Verdana"
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

function bar_plot(metab,df)
    # Handle some special cases
    if metab == "H_x"
        fmt = ".1e"
    elseif metab == "K_x"
        fmt = ".0f"
    elseif metab in ["Cr","Cr_i","ATP_e","ATP_i"]
        fmt = ".1f"
    else
        fmt = ".2"
    end

    df_plot = filter(:metabolite => ==(metab),df)
    return df_plot |> @vlplot(
        title={
            text=metab,
            fontSize=7,
            fontWeight=:normal,
            fontStyle=:italic,
            offset=2
        },
        x={
            "concentration:q", 
            axis={
                title="Conc. (mM)", 
                titleFontWeight=:normal, 
                titleFontSize=6, 
                grid=false,
                tickSize=0,
                labels=false,
                domainColor=:black,
                titlePadding=0
            },
        scale={domain=[0,1.4*maximum(df_plot.concentration)]},
        },
        y={
            "group:n", 
            axis={
                title="",
                labels=false,
                grid=false,
                tickSize=0,
                domainColor=:black
            }, 
            sort=:descending
        },
        config={
            view={stroke=:transparent},
            legend={title=""},
            font="Verdana"
        }
    ) +
    @vlplot(
        mark={:bar, height=5.67}, width=40, height=18,
        color={"group:n", scale={range=["#009DA3","#F98F00"]}, sort=:descending, legend=nothing}
    ) +
    @vlplot(
        mark={
            :text,
            align=:left,
            baseline=:middle,
            fontSize=6,
            dx=1
        },
        text={"concentration:q",format=fmt}
    )
end


## 
md"""
Generate results at low workloads
"""
VO2_human = 6.05e-4
df_results = simulate_metabolomics(VO2_human)
CSV.write("../output/oxphos_model_fc_low_wl.csv", df_results)
fig = plot_metabolomics_results(df_results)
savevl(fig,"../output/metabolites_low_wl")
enlarge(fig,2)

## 
md"""
Generate results at high workloads
"""
df_results = simulate_metabolomics(2*VO2_human)
CSV.write("../output/oxphos_model_fc_high_wl.csv", df_results)
fig = plot_metabolomics_results(df_results)
savevl(fig,"../output/metabolites_high_wl")
enlarge(fig,2)

##
md"""
Run a sensitivity analysis of the figures with respect to the parameters μ_r, μ_DH, and μa_ANT
"""
sys = ageing_sys()

function perturbed_simulations(sys,param,ratio)
    param_store = sys.defaults[param]
    sys.defaults[param] = sys.defaults[param] + L(ratio)
    atpase_params = workload_param_vec(max_wl=10)
    array_sols = [
        :Young => simulate_FCs(sys,0,p=atpase_params), 
        :Old => simulate_FCs(sys,1,p=atpase_params)
    ];
    sys.defaults[param] = param_store

    dfs = [extract_vars(sols,string(scale)) for (scale,sols) in array_sols]
    df = vcat(dfs...)
    return df
end

@parameters μ_r μa_DH μa_ANT
function simulate_and_plot(sys,param,ratio)
    df = perturbed_simulations(sys,param,ratio)
    fig_PCrATP = PCrATP_plot(df)
    fig_ADP = ADP_plot(df)
    return (df,fig_PCrATP,fig_ADP)
end

(df_μ_r,fig_μ_r_PCrATP,fig_μ_r_ADP) = simulate_and_plot(sys,μ_r,0.7)
(df_μa_DH,fig_μa_DH_PCrATP,fig_μa_DH_ADP) = simulate_and_plot(sys,μa_DH,0.7)
(df_μa_ANT,fig_μa_ANT_PCrATP,fig_μa_ANT_ADP) = simulate_and_plot(sys,μa_ANT,0.7)
savevl(fig_μ_r_PCrATP,"../output/sensitivity/PCrATP_mu_r_down")
savevl(fig_μ_r_ADP,"../output/sensitivity/ADP_mu_r_down")
savevl(fig_μa_DH_PCrATP,"../output/sensitivity/PCrATP_DH_down")
savevl(fig_μa_DH_ADP,"../output/sensitivity/ADP_DH_down")
savevl(fig_μa_ANT_PCrATP,"../output/sensitivity/PCrATP_ANT_down")
savevl(fig_μa_ANT_ADP,"../output/sensitivity/ADP_ANT_down")

(df_μ_r2,fig_μ_r_PCrATP2,fig_μ_r_ADP2) = simulate_and_plot(sys,μ_r,1.3)
(df_μa_DH2,fig_μa_DH_PCrATP2,fig_μa_DH_ADP2) = simulate_and_plot(sys,μa_DH,1.3)
(df_μa_ANT2,fig_μa_ANT_PCrATP2,fig_μa_ANT_ADP2) = simulate_and_plot(sys,μa_ANT,1.3)
savevl(fig_μ_r_PCrATP2,"../output/sensitivity/PCrATP_mu_r_up")
savevl(fig_μ_r_ADP2,"../output/sensitivity/ADP_mu_r_up")
savevl(fig_μa_DH_PCrATP2,"../output/sensitivity/PCrATP_DH_up")
savevl(fig_μa_DH_ADP2,"../output/sensitivity/ADP_DH_up")
savevl(fig_μa_ANT_PCrATP2,"../output/sensitivity/PCrATP_ANT_up")
savevl(fig_μa_ANT_ADP2,"../output/sensitivity/ADP_ANT_up")
