using ModelingToolkit, DifferentialEquations, Sundials, .OxPhos

@parameters t
@parameters μa_ATPase μa_ANT

# Import the model for use in ageing heart
function ageing_sys()
    @variables NAD_x(t) NADH_x(t)
    sys = bg_model_invivo(stat=[NAD_x,NADH_x]) # Set NAD and NADH to be constant concentrations
    sys.defaults[μa_ANT] = OxPhos.updated_μa_ANT() # Use fitted parameter for ANT
    return sys
end

# Set model options
tspan = (0.0,1000.0)
options = Dict(:reltol => 1e-7, :abstol => 1e-10, :alg => CVODE_BDF())

# Define parameter range for ATP consumption
function ATPase_params(;n_steps=100,stop=7)
    ATP_consumption_ratio = range(0.1,stop=stop,length=n_steps)
    ps = OxPhos.parameters_bg_invivo()
    return @. ps[μa_ATPase] + L(ATP_consumption_ratio)
end

@parameters μa_ANT
workload_param_vec(;n_steps=100,max_wl=7) = [[μa_ATPase => m] for m in ATPase_params(n_steps=n_steps,stop=max_wl)]

# Define function to apply fold changes in the data
protein_FC(r) = read_entry(df_averaged_FC,:name,r).FC[1]
adjust_logFC(μ,FC,scale=1) = μ + L(1 + scale*(FC-1))

function generate_p_FC(scale=1)
    @parameters μa_C1 μa_C3 μa_C4 μa_F1 μa_ANT μa_KH μa_PiHt
    p0 = OxPhos.parameters_bg_invivo()
    return Dict(
        μa_C1 => adjust_logFC(p0[μa_C1],protein_FC("Mitochondrial complex I"),scale),
        μa_C3 => adjust_logFC(p0[μa_C3],protein_FC("Mitochondrial complex III"),scale),
        μa_C4 => adjust_logFC(p0[μa_C4],protein_FC("Mitochondrial complex IV"),scale),
        μa_F1 => adjust_logFC(p0[μa_F1],protein_FC("Mitochondrial complex V"),scale),
        μa_ANT => adjust_logFC(p0[μa_ANT],protein_FC("Adenine nucleotide translocator"),scale),
        μa_KH => adjust_logFC(p0[μa_KH],protein_FC("K/H antiport"),scale),
        μa_PiHt => adjust_logFC(p0[μa_PiHt],protein_FC("Phosphate carrier protein"),scale)
    )
end
metabolite_FC(m) = read_entry(df_metabolite_FC,:name,m).FC[1]
read_entry(df,field,value) = df[getproperty(df,field) .== value, :]

adjust_FC(u,FC,scale=1) = u * (1 + scale*(FC-1))
function adjust_u0(scale=1,Cr_scale=0) 
    u0 = OxPhos.u0_bg_invivo()
    @variables Cr(t) PCr(t) Cr_i(t) PCr_i(t) NAD_x(t) NADH_x(t)
    return [
        NAD_x => adjust_FC(u0[NAD_x],metabolite_FC("Nicotinamide adenine dinucleotide"),scale),
        NADH_x => adjust_FC(u0[NADH_x],metabolite_FC("Reduced nicotinamide Adenine Dinucleotide"),scale),
        # Assume that creatine fold change also applies to phosphocreatine
        Cr => adjust_FC(u0[Cr],metabolite_FC("Creatine"),Cr_scale),
        PCr => adjust_FC(u0[PCr],metabolite_FC("Phosphocreatine"),Cr_scale),
        Cr_i => adjust_FC(u0[Cr_i],metabolite_FC("Creatine"),Cr_scale),
        PCr_i => adjust_FC(u0[PCr_i],metabolite_FC("Phosphocreatine"),Cr_scale),
    ]
end

function simulate_FCs(sys,scale=1,Cr_scale=0; p=workload_param_vec())
    u0 = adjust_u0(scale,Cr_scale)
    sols = multi_sim(sys,tspan,p;u0=u0,options...)
end
