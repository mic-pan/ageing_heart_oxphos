module OxPhos

using ModelingToolkit, Symbolics, Optimization, OptimizationOptimJL
using JSON, Sundials, Interpolations
export val, E, L, μ
export beard_model, bg_model_invitro, bg_model_invivo
include("simulation_functions.jl")

@parameters t
D = Differential(t)
RT_default = 2.4734

# Useful functions
val(s) = s.val.metadata[Symbolics.VariableDefaultValue]
E(μ,RT=RT_default) = exp(μ/RT)
L(x,RT=RT_default) = RT*log(x)
μ(μ0,c,RT=RT_default) = μ0 + L(c,RT)

# Define parameters and variables for Beard model
@parameters Qtot dG_C3o k_Pi3 k_Pi4 x_C1 x_C3 gamma p_A AMP_e x_MgA K_DT x_DH r k_Pi1 k_Pi2 n_A C_im NADtot F RT ATP_e x_buff x_ANT k_mADP Ctot dG_C1o x_C x_AK K_AK ADP_e K_DD Mg_tot x_KH x_Pi1 k_PiH dG_C4o x_C4 k_O2 pH_e W_m K_e x_F1 dG_F1o x_Hle Pi_e x_Pi2 x_K k_dHPi k_dHatp k_dHadp
@variables Cox(t) dG_C3op(t) dPsi(t) H_x(t) J_C3(t) Pi_x(t) Q(t) QH2(t) dG_H(t) Cred(t) AMP_i(t) J_AKi(t) J_AMP(t) W_i(t) ADP_i(t) J_ADP(t) J_AKi(t) J_ANT(t) ATP_fx(t) ATP_x(t) ATP_mx(t) J_MgATPx(t) Mg_x(t) J_DH(t) NAD_x(t) NADH_x(t) J_C1(t) J_C4(t) J_ANT(t) J_Hle(t) J_K(t) J_F1(t) H_i(t) W_x(t) J_ATP(t) ATP_i(t) J_Pi1(t) J_DH(t) J_KH(t) Psi_x(t) Psi_i(t) ADP_fi(t) ATP_fi(t) ADP_fx(t) ATP_mi(t) J_MgATPi(t) Mg_i(t) dG_C1op(t) ADP_fe(t) ADP_me(t) K_i(t) K_x(t) O2(t) Mg_e(t) Mg_i(t) ADP_mi(t) J_MgADPi(t) ADP_mx(t) J_MgADPx(t) H2PIi(t) Pi_i(t) H2PIx(t) dG_C4op(t) H_e(t) ADP_x(t) J_Pi2(t)

"""
Return a dictionary with the parameters of the Beard model.
Note: x_AK set to 1e6 so that reaction is at equilibrium
"""
parameters_beard() = Dict(
    Qtot => 0.00135, 
    dG_C3o => -32.53, 
    k_Pi3 => 0.00019172, 
    k_Pi4 => 0.02531, 
    x_C3 => 0.091737, 
    gamma => 5.99, 
    p_A => 85.0, 
    AMP_e => 0.0, 
    x_MgA => 1.0e6, 
    K_DT => 2.4e-5, 
    x_DH => 0.09183, 
    r => 4.5807, 
    k_Pi1 => 0.00013413, 
    k_Pi2 => 0.00067668, 
    n_A => 3, 
    C_im => 6.756756756756757e-6, 
    NADtot => 0.00297, 
    F => 0.096484, 
    RT => RT_default, 
    ATP_e => 0.0, 
    x_buff => 100.0, 
    x_ANT => 0.0079204, 
    k_mADP => 3.5e-6, 
    Ctot => 0.0027,
    dG_C1o => -69.37, 
    x_C1 => 0.36923, 
    x_AK => 1e6,
    K_AK => 0.4331, 
    ADP_e => 0.0013, 
    K_DD => 0.000347, 
    Mg_tot => 0.005, 
    x_KH => 2.9802e7, 
    x_Pi1 => 339430.0, 
    k_PiH => 0.00045082, 
    dG_C4o => -122.94, 
    x_C4 => 3.2562e-5, 
    k_O2 => 0.00012, 
    pH_e => 7.1, 
    W_m => 0.72376, 
    K_e => 0.15, 
    x_F1 => 150.93, 
    dG_F1o => 36.03, 
    x_Hle => 250.0, 
    Pi_e => 0.005, 
    x_Pi2 => 327.0, 
    x_K => 0.0, 
    k_dHPi => 1.7782794100389227e-7, 
    k_dHatp => 3.311311214825908e-7, 
    k_dHadp => 5.128613839913648e-7
)


"""
Return a dictionary with the initial conditions of the Beard model.
Small concentrations are used for ADP_i, ATP_i, ATP_mi and ADP_mi to avoid taking logs of zero. 
"""
u0_beard() = Dict(
    dPsi => 160.0, 
    H_x => 6.30957344480193e-8, 
    Pi_x => 0.001, 
    QH2 => 0.0008, 
    Cred => 0.001, 
    AMP_i => 0.0, 
    ADP_i => 2e-15, 
    ATP_x => 0.0, 
    ATP_mx => 0.0, 
    Mg_x => 0.005, 
    NADH_x => 0.0015, 
    ATP_i => 2e-15, 
    ATP_mi => 1e-15, 
    K_x => 0.14, 
    O2 => 2.6e-5, 
    ADP_mx => 0.0, 
    ADP_mi => 1e-15, 
    Pi_i => 0.001, 
    ADP_x => 0.01
)

"""
Return an ODE system of the Beard model.
Ref: https://doi.org/10.1371/journal.pcbi.0010036
"""
function beard_model()
    eqns_beard = Equation[
        Q ~ Qtot - QH2,
        dG_C3op ~ dG_C3o + 2.0RT*log(1.0e7H_x) - RT*log(QH2 / Q),
        J_C3 ~ (x_C3*(1.0 + Pi_x / k_Pi3)*(Cox*exp((2.0F*dPsi - dG_C3op - 4.0dG_H) / (2.0RT)) - Cred)) / (1.0 + Pi_x / k_Pi4),
        D(AMP_i) ~ (J_AKi + J_AMP) / W_i,
        J_AMP ~ gamma*p_A*(AMP_e - AMP_i),
        D(ADP_i) ~ (J_ADP - 2.0J_AKi - J_ANT) / W_i,
        ATP_fx ~ ATP_x - ATP_mx,
        J_MgATPx ~ x_MgA*(ATP_fx*Mg_x - K_DT*ATP_mx),
        J_DH ~ (x_DH*(r*NAD_x - NADH_x)*(1.0 + Pi_x / k_Pi1)) / (1.0 + Pi_x / k_Pi2),
        D(dPsi) ~ (2.0J_C3 + 4.0J_C1 + 4.0J_C4 - J_ANT - J_Hle - J_K - n_A*J_F1) / C_im,
        NAD_x ~ NADtot - NADH_x,
        J_ADP ~ gamma*p_A*(ADP_e - ADP_i),
        D(ATP_i) ~ (J_AKi + J_ANT + J_ATP) / W_i,
        J_Hle ~ (x_Hle*(H_i*exp((F*dPsi) / RT) - H_x)*dPsi) / (exp((F*dPsi) / RT) - 1.0),
        J_Pi2 ~ gamma*x_Pi2*(Pi_e - Pi_i),
        D(Cred) ~ (2.0J_C3 - 2.0J_C4) / W_i,
        D(Pi_i) ~ (J_Pi2 - J_Pi1) / W_i,
        J_K ~ (x_K*(K_i*exp((F*dPsi) / RT) - K_x)*dPsi) / (exp((F*dPsi) / RT) - 1.0),
        dG_H ~ F*dPsi + RT*log(H_i / H_x),
        D(NADH_x) ~ (J_DH - J_C1) / W_x,
        J_ATP ~ gamma*p_A*(ATP_e - ATP_i),
        D(H_x) ~ (x_buff*((n_A - 1.0)*J_F1 + 2.0J_Pi1 + J_DH + J_Hle - 5.0J_C1 - J_KH - 2.0J_C3 - 4.0J_C4)*H_x) / W_x,
        Psi_x ~ -0.65dPsi,
        Psi_i ~ 0.35dPsi,
        J_ANT ~ (x_ANT*((-ADP_fx) / (ATP_fx*exp((-F*Psi_x) / RT) + ADP_fx) + ADP_fi / (ATP_fi*exp((-F*Psi_i) / RT) + ADP_fi))*ADP_fi) / (k_mADP + ADP_fi),
        D(ATP_mi) ~ J_MgATPi / W_i,
        ATP_fi ~ ATP_i - ATP_mi,
        J_MgATPi ~ x_MgA*(ATP_fi*Mg_i - K_DT*ATP_mi),
        Cox ~ Ctot - Cred,
        dG_C1op ~ dG_C1o - RT*log(1.0e7H_x) - RT*log(Q / QH2),
        J_C1 ~ x_C1*(NADH_x*exp((-dG_C1op - 4.0dG_H) / RT) - NAD_x),
        J_AKi ~ x_AK*(K_AK*(ADP_i^2) - AMP_i*ATP_i),
        ADP_fe ~ ADP_e - ADP_me,
        ADP_me ~ 0.5ADP_e + 0.5K_DD + 0.5Mg_tot - 0.5sqrt((ADP_e + K_DD + Mg_tot)^2.0 - 4.0ADP_e*Mg_tot),
        D(ATP_x) ~ (J_F1 - J_ANT) / W_x,
        J_KH ~ x_KH*(H_x*K_i - H_i*K_x),
        D(Pi_x) ~ (J_Pi1 - J_F1) / W_x,
        D(O2) ~ 0.0,
        D(Mg_x) ~ (-J_MgADPx - J_MgATPx) / W_x,
        Mg_e ~ Mg_tot - ADP_me,
        Mg_i ~ Mg_e,
        D(ATP_mx) ~ J_MgATPx / W_x,
        ADP_fi ~ ADP_i - ADP_mi,
        J_MgADPi ~ x_MgA*(ADP_fi*Mg_i - K_DD*ADP_mi),
        D(K_x) ~ (J_K + J_KH) / W_x,
        D(ADP_mx) ~ J_MgADPx / W_x,
        H2PIi ~ (H_i*Pi_i) / (H_i + k_dHPi),
        H2PIx ~ (H_x*Pi_x) / (H_x + k_dHPi),
        J_Pi1 ~ (x_Pi1*(H2PIi*H_x - H2PIx*H_i)) / (k_PiH + H2PIi),
        dG_C4op ~ dG_C4o - 0.5RT*log(O2) - 2.0RT*log(1.0e7H_x),
        J_C4 ~ (x_C4*(Cred*exp((-dG_C4op - 2.0dG_H) / (2.0RT)) - Cox*exp((F*dPsi) / RT))*Cred) / (Ctot*(1.0 + k_O2 / O2)),
        D(QH2) ~ (J_C1 - J_C3) / W_x,
        H_e ~ 10.0^(-pH_e),
        W_x ~ 0.9W_m,
        W_i ~ 0.1W_m,
        H_i ~ H_e,
        K_i ~ K_e,
        D(ADP_x) ~ (J_ANT - J_F1) / W_x,
        J_F1 ~ x_F1*((K_DD*ADP_mx*Pi_x*exp((n_A*dG_H - dG_F1o) / RT)) / K_DT - ATP_mx),
        ADP_fx ~ ADP_x - ADP_mx,
        J_MgADPx ~ x_MgA*(ADP_fx*Mg_x - K_DD*ADP_mx),
        D(ADP_mi) ~ J_MgADPi / W_i
    ]

    p_vals = parameters_beard()
    ps = collect(keys(p_vals))
    sts_beard = [
        dPsi, H_x, Pi_x, QH2, Cred, AMP_i, ADP_i, ATP_x, ATP_mx, Mg_x, NADH_x, ATP_i, ATP_mi, K_x, O2, ADP_mi, ADP_mx, Pi_i, ADP_x
    ]

    odesys = ODESystem(eqns_beard,t,sts_beard,ps;name=:oxphos,defaults=merge(p_vals,u0_beard()))
    sys = structural_simplify(odesys)
    return sys
end

# Define the parameters and variables used in the bond graph model
@parameters μ0_NADH_x μ0_NAD_x μ0_Q μ0_QH2 μ0_H μ0_Cox μ0_Cred μ_H2O μ0_O2 μ0_fATP μ0_fADP μ0_Pi μ0_AMP μ0_Mg μ0_mATP μ0_mADP μ0_K μ0_B μ0_HB
@parameters μa_C1 μa_C3 μa_C4 μa_F1 μa_MgATP μa_MgADP μa_Hle μa_K μa_KH μa_DH μ_r μa_fATPt μa_fADPt μa_mATPt μa_mADPt μa_AMPt μa_Pit μa_AK μa_ANT μa_PiHt μa_Hb
@parameters μreg_Pi1 μreg_Pi2 μreg_Pi3 μreg_Pi4 μreg_O2 μreg_Cred μd_HPi μreg_PiHt μreg_ANT
@parameters Btot K_dHbuff
@variables μ_NADH_x(t) μ_Q(t) μ_H_i(t) μ_H_x(t) μ_NAD_x(t) μ_QH2(t) μ_Cox(t) μ_Cred(t) μ_O2(t) μ_ATP_fx(t) μ_ATP_mx(t) μ_ATP_fi(t) μ_ATP_mi(t) μ_ADP_fx(t) μ_ADP_mx(t) μ_ADP_fi(t) μ_ADP_mi(t) μ_Pi_x(t) μ_Pi_i(t) μ_Mg_x(t) μ_Mg_i(t) μ_K_i(t) μ_K_x(t) μ_ATP_fe(t) μ_ATP_me(t) μ_ADP_fe(t) μ_ADP_me(t) μ_AMP_e(t) μ_AMP_i(t) μ_Pi_e(t) ATP_fe(t) ATP_me(t) J_fATP(t) J_mATP(t) J_fADP(t) J_mADP(t) J_Hbuff(t) B_x(t) HB_x(t) μ_B_x(t) μ_HB_x(t)

"""
Returns a dict with the bond graph species parameters (corresponding to standard free energies of formation).

Unless otherwise stated, parameters were taken from Wu et al. (2007; https://doi.org/10.1074/jbc.M701024200). 
μ0_Mg was taken from eQuilibrator with ionic strength 0.17.
μ0_K was assumed to be zero.
Because the buffer is a phenomenological species, it is assigned a standard free energy of zero. The parameter for the H-bound buffer is calculated using the dissociation constant.
"""
function bg_species_parameters()
    v_μ0_Mg = -458.2
    v_μ0_fATP = -2771.00
    v_μ0_fADP = -1903.96
    v_μ0_H = 0.0
    v_μ0_B = 0.0

    pb = parameters_beard()
    pbuff = buffering_parameters()

    return Dict(
        μ0_NADH_x => 39.31, 
        μ0_NAD_x => 18.10, 
        μ0_Q => 65.17, 
        μ0_QH2 => -23.30, 
        μ0_H => v_μ0_H, 
        μ0_Cox => -6.52, 
        μ0_Cred => -27.41, 
        μ_H2O => -235.74, 
        μ0_O2 => 16.40, 
        μ0_fATP => v_μ0_fATP, 
        μ0_fADP => v_μ0_fADP, 
        μ0_Pi => -1098.27, 
        μ0_AMP => -1034.66,
        μ0_Mg => v_μ0_Mg,
        μ0_mATP => v_μ0_fATP+v_μ0_Mg+L(pb[K_DT]),
        μ0_mADP => v_μ0_fADP+v_μ0_Mg+L(pb[K_DD]),
        μ0_K => 0.0,
        μ0_B => v_μ0_B, 
        μ0_HB => v_μ0_H + v_μ0_B + L(pbuff[K_dHbuff])
    )
end

"""
Returns a dict with the reaction parameters (roughly corresponding the the activation energy of reaction).

In most cases, the reaction parameter was chosen to match the forward rate constant of the kinetic model.
The buffering reaction was assumed to be rapid compared to the other reactions in the system.
"""
function bg_reaction_parameters()
    p_be = parameters_beard()
    p_μ0 = bg_species_parameters()

    H_nom = 10^(-7.2)
    μ_H_nom = p_μ0[μ0_H] + L(H_nom)
    μa_DH_val = L(p_be[x_DH])-p_μ0[μ0_NADH_x]-μ_H_nom
    μ_r_val = L(p_be[x_DH]*p_be[r])-μa_DH_val-p_μ0[μ0_NAD_x]

    return Dict(
        μa_C1 => L(1e7*p_be[x_C1])-p_be[dG_C1o]-p_μ0[μ0_H]+p_μ0[μ0_QH2]-p_μ0[μ0_Q]-p_μ0[μ0_NADH_x],
        μa_C3 => L(1e-7*p_be[x_C3])+(-p_be[dG_C3o]+p_μ0[μ0_Q]-p_μ0[μ0_QH2]-2p_μ0[μ0_Cox]+2p_μ0[μ0_H])/2,
        μa_C4 => L(1e7*p_be[x_C4])-(p_be[dG_C4o]+2p_μ0[μ0_Cred]+1/2*p_μ0[μ0_O2])/2,
        μa_F1 => L(p_be[x_F1]*p_be[K_DD]/p_be[K_DT]/H_nom)-p_be[dG_F1o]-p_μ0[μ0_mADP]-p_μ0[μ0_Pi]-p_μ0[μ0_H],
        μa_MgATP => L(p_be[x_MgA])-p_μ0[μ0_Mg]-p_μ0[μ0_fATP],
        μa_MgADP => L(p_be[x_MgA])-p_μ0[μ0_Mg]-p_μ0[μ0_fADP],
        μa_Hle => L(p_be[x_Hle])-p_μ0[μ0_H],
        #μa_K => L(p_be[x_K])-p_μ0[μ0_K],
        μa_KH => L(p_be[x_KH])-p_μ0[μ0_H]-p_μ0[μ0_K],
        μa_DH => μa_DH_val,
        μ_r => μ_r_val,
        μa_fATPt => L(p_be[gamma]*p_be[p_A])-p_μ0[μ0_fATP],
        μa_fADPt => L(p_be[gamma]*p_be[p_A])-p_μ0[μ0_fADP],
        μa_mATPt => L(p_be[gamma]*p_be[p_A])-p_μ0[μ0_mATP],
        μa_mADPt => L(p_be[gamma]*p_be[p_A])-p_μ0[μ0_mADP],
        μa_AMPt => L(p_be[gamma]*p_be[p_A])-p_μ0[μ0_AMP],
        μa_Pit => L(p_be[gamma]*p_be[x_Pi2])-p_μ0[μ0_Pi],
        μa_AK => L(p_be[x_AK])-2p_μ0[μ0_fADP],
        μa_ANT => L(p_be[x_ANT]),
        μa_PiHt => L(p_be[x_Pi1]/p_be[k_dHPi]/p_be[k_PiH])-p_μ0[μ0_Pi]-2p_μ0[μ0_H],
        μa_Hb => L(1e9) - p_μ0[μ0_H] - p_μ0[μ0_B]
    )
end

"""
Return a dict with the phenomenological regulation parameters. The parameters were chosen to match the equations in the Beard model.
"""
function bg_regulation_parameters()
    p_be = parameters_beard()
    p_μ0 = bg_species_parameters()

    return Dict(
        μreg_Pi1 => p_μ0[μ0_Pi] + L(p_be[k_Pi1]),
        μreg_Pi2 => p_μ0[μ0_Pi] + L(p_be[k_Pi2]),
        μreg_Pi3 => p_μ0[μ0_Pi] + L(p_be[k_Pi3]),
        μreg_Pi4 => p_μ0[μ0_Pi] + L(p_be[k_Pi4]),
        μreg_O2 => p_μ0[μ0_O2] + L(p_be[k_O2]),
        μreg_Cred => p_μ0[μ0_Cred] + L(p_be[Ctot]),
        μd_HPi => p_μ0[μ0_H] + L(p_be[k_dHPi]),
        μreg_PiHt => p_μ0[μ0_Pi] + L(p_be[k_PiH]),
        μreg_ANT => p_μ0[μ0_fADP] + L(p_be[k_mADP])
    )
end

"""
Optimise for the buffering parameters.

The Beard model uses the term x_buff*H in dH_x/dt, which is not biochemically plausible. Instead, the term 1/(1 + Bt*K/(H+K)^2) was used instead, which is derived from H + B ⇄ HB at equilibrium. The updated equations ensures that total H is conserved.

The parameters of the new equations were chosen to give the best fit to the Beard equation from pH 7.0 to 7.4. A sum of squares cost function was used, with a regularisation term.
"""
function fit_buffering_parameters()
    @variables Bt K
    x_buff = 100

    f1(H) = x_buff*H # function in Beard model
    f2(H) = 1/(1 + Bt*K/(H+K)^2) # equation derived from the biophysics of H + B ⇄ HB at equilibrium
    diff(H) = f1(H) - f2(H)

    H_fit = 10 .^ (-7.4:0.01:-7.0)
    loss = sum(diff(H)^2 for H in H_fit)

    sys = OptimizationSystem(loss,[Bt,K],[];name=:buffer)
    u0 = [Bt => 1.0, K => 1.0]
    prob = OptimizationProblem(sys,u0,[],grad=true,hess=true)
    u_opt = solve(prob,NewtonTrustRegion())
    return u_opt
end

"""
Return the parameters associated with buffering. By default, the the results of previously run optimisation code is returned. The run_opt option can be set to true to run the optimisation again.
"""
function buffering_parameters(;run_opt=false)
    if run_opt
        (Bt_opt, K_opt) = fit_buffering_parameters()
        return Dict(
            Btot => Bt_opt, 
            K_dHbuff => K_opt
        )
    else
        return Dict(
            Btot => 0.04071009151505011, 
            K_dHbuff => 6.689164053279534e-8
        )
    end
end

@variables ΔG_C1(t) ΔG_C3(t) ΔG_C4(t) ΔG_F1(t) ΔG_Pi1(t)

"""
The volume ratio between the intermembrane space and mitochondria.
"""
ims_mito_ratio() = 0.1

"""
The volume ratio between the mitochondrial matrix and entire mitochondria.
"""
matrix_mito_ratio() = 1-ims_mito_ratio()

"""
Returns the cytosolic/mitochondrial volume ratio for the given species
"""
function cyto_mito_ratio(species=:human)
    if species == :rat
        return 0.399/0.143 # From Table 5 of Vinnakota and Bassingthwaighte (2004; https://doi.org/10.1152/ajpheart.00478.2003). The scaling constant is the ratio of cytosolic to mitochondrial volume.
    elseif species == :human
        volfrac_mito = 0.2533
        return (1-volfrac_mito)/volfrac_mito # Table 1, Barth et al. (1992; https://doi.org/10.1016/0022-2828(92)93381-S) 
    end
    return 
end

"""
Equations shared between all versions of the bond graph model.

Notes on J_Pi1:
1. The protons in the original PiHt flux J_Pi1 = (x_Pi1*(H2PIi*H_x - H2PIx*H_i)) / (k_PiH + H2PIi) were assigned incorrectly as the transporter is a co-transporter. This has been fixed here.
2. Since H2PO4 wasn't a state variable in the original model, and J_Pi1 contributes to the rate of Pi rather than H2PO4, I've modified J_Pi1 to use Pi as a species. So the reaction is now Pi_i + 2H_i ⇌ Pi_x + 2H_x.
3. Since the regulation term 1/(H_i + k_dHPi) isn't dimensionless, I've equated it to (1/k_dHPi)*1/(1+H_i/k_dHPi), where the first part of the product is absorbed into the rate constant μa_PiHt. The same applies for k_PiH.
"""
equations_bg_base() = Equation[
    # Ce:NADH_x
    D(NADH_x) ~ (J_DH - J_C1) / W_x,
    μ_NADH_x ~ μ(μ0_NADH_x,NADH_x),

    # Ce:NAD_x
    #NAD_x ~ NADtot - NADH_x, 
    D(NAD_x) ~ -(J_DH - J_C1) / W_x,
    μ_NAD_x ~ μ(μ0_NAD_x,NAD_x),

    # Ce:Q
    Q ~ Qtot - QH2, # D(Q) ~ -(J_C1 - J_C3) / W_x
    μ_Q ~ μ(μ0_Q,Q),

    # Ce:QH2
    D(QH2) ~ (J_C1 - J_C3) / W_x,
    μ_QH2 ~ μ(μ0_QH2,QH2),

    # Ce:H_x
    D(H_x) ~ ((n_A-1)*J_F1 + 2J_Pi1 + J_DH + J_Hle - 5J_C1 - J_KH - 2J_C3 - 4J_C4 - J_Hbuff) / W_x,
    μ_H_x ~ μ(μ0_H,H_x),
    
    #Ce:B_x
    D(B_x) ~ -J_Hbuff/W_x,
    μ_B_x ~ μ(μ0_B,B_x),

    #Ce:HB_x
    D(HB_x) ~ J_Hbuff/W_x,
    μ_HB_x ~ μ(μ0_HB,HB_x),

    # Se:H_i
    H_e ~ 10.0^(-pH_e),
    H_i ~ H_e,
    μ_H_i ~ μ(μ0_H,H_i),

    # Ce:Cox
    Cox ~ Ctot - Cred, # D(Cox) ~ -(2J_C3 - 2J_C4) / W_i
    μ_Cox ~ μ(μ0_Cox,Cox),

    # Ce:Cred
    D(Cred) ~ (2J_C3 - 2J_C4) / W_i,
    μ_Cred ~ μ(μ0_Cred,Cred),

    # Se:O2
    D(O2) ~ 0.0, # The concentration of oxygen is constant
    μ_O2 ~ μ(μ0_O2,O2),

    # Ce:ATP_fx
    D(ATP_fx) ~ (-J_ANT - J_MgATPx) / W_x, #ATP_fx ~ ATP_x - ATP_mx,
    μ_ATP_fx ~ μ(μ0_fATP,ATP_fx),

    # Ce:ATP_mx
    D(ATP_mx) ~ (J_F1 + J_MgATPx) / W_x,
    μ_ATP_mx ~ μ(μ0_mATP,ATP_mx),

    # Ce:ATP_fi
    D(ATP_fi) ~ (J_AKi + J_ANT + J_fATP - J_MgATPi) / W_i, # ATP_fi ~ ATP_i - ATP_mi,
    μ_ATP_fi ~ μ(μ0_fATP,ATP_fi), 

    # Ce:ADP_fx
    D(ADP_fx) ~ (J_ANT - J_MgADPx) / W_x, #ADP_fx ~ ADP_x - ADP_mx,
    μ_ADP_fx ~ μ(μ0_fADP,ADP_fx),

    # Ce:ADP_mx
    D(ADP_mx) ~ (J_MgADPx - J_F1) / W_x,
    μ_ADP_mx ~ μ(μ0_mADP,ADP_mx),

    # Ce:ADP_fi
    D(ADP_fi) ~ (J_fADP - 2J_AKi - J_ANT - J_MgADPi) / W_i, # ADP_fi ~ ADP_i - ADP_mi,
    μ_ADP_fi ~ μ(μ0_fADP,ADP_fi),

    # Ce:Pi_x
    D(Pi_x) ~ (J_Pi1 - J_F1) / W_x,
    μ_Pi_x ~ μ(μ0_Pi,Pi_x),

    # Ce:Pi_i
    D(Pi_i) ~ (J_Pi2 - J_Pi1) / W_i,
    μ_Pi_i ~ μ(μ0_Pi,Pi_i),

    # Ce:Mg_x
    D(Mg_x) ~ (-J_MgADPx - J_MgATPx) / W_x,
    μ_Mg_x ~ μ(μ0_Mg,Mg_x),

    # Se:K_i 
    K_i ~ K_e,
    μ_K_i ~ μ(μ0_K,K_i),

    # Ce:K_x
    D(K_x) ~ (J_K + J_KH) / W_x,
    μ_K_x ~ μ(μ0_K,K_x),

    # Ce:AMP_i
    D(AMP_i) ~ (J_AKi + J_AMP) / W_i,
    μ_AMP_i ~ μ(μ0_AMP,AMP_i),

    # R:Hbuff
    J_Hbuff ~ E(μa_Hb + μ_H_x + μ_B_x) - E(μa_Hb + μ_HB_x),

    # R:C1
    J_C1 ~ E(μa_C1 + μ_NADH_x + μ_Q - μ_QH2 + μ_H_x - 4dG_H) - E(μa_C1 + μ_NAD_x),
    
    # R:C3
    J_C3 ~ (E(μa_C3 + (2μ_Cox - 4dG_H - 2μ_H_x + 2F*dPsi + μ_QH2 - μ_Q)/2) - E(μa_C3 + μ_Cred))*(1 + E(μ_Pi_x - μreg_Pi3))/ (1 + E(μ_Pi_x - μreg_Pi4)),

    # R:C4
    J_C4 ~ (E(μa_C4 + (2μ_H_x + 2μ_Cred + 1/2*μ_O2 - 2dG_H)/2) - E(μa_C4 + (2F*dPsi + 2μ_Cox + μ_H2O)/2))*1/(1 + E(μreg_O2 - μ_O2))*E(μ_Cred-μreg_Cred),

    # R:F1
    J_F1 ~ E(μa_F1 + μ_ADP_mx + μ_Pi_x + μ_H_x + n_A*dG_H) - E(μa_F1 + μ_ATP_mx + μ_H2O),

    # R:J_MgATPx
    J_MgATPx ~ E(μa_MgATP + μ_Mg_x + μ_ATP_fx) - E(μa_MgATP + μ_ATP_mx),
    # R:J_MgATPi
    J_MgATPi ~ E(μa_MgATP + μ_Mg_i + μ_ATP_fi) - E(μa_MgATP + μ_ATP_mi),
    # R:J_MgADPx
    J_MgADPx ~ E(μa_MgADP + μ_Mg_x + μ_ADP_fx) - E(μa_MgADP + μ_ADP_mx),
    # R:J_MgADPi
    J_MgADPi ~ E(μa_MgADP + μ_Mg_i + μ_ADP_fi) - E(μa_MgADP + μ_ADP_mi),

    # R:Hle
    J_Hle ~ (E(μa_Hle+μ_H_i+F*dPsi) - E(μa_Hle+μ_H_x))*dPsi/(exp((F*dPsi) / RT) - 1),

    # R:J_K
    J_K ~ 0, # Set to zero since x_K is zero

    # R:KH
    J_KH ~ E(μa_KH + μ_H_x + μ_K_i) - E(μa_KH + μ_H_i + μ_K_x),

    # R: DH
    J_DH ~ (E(μa_DH+μ_r+μ_NAD_x) - E(μa_DH+μ_NADH_x+μ_H_x))*(1 + E(μ_Pi_x - μreg_Pi1))/(1 + E(μ_Pi_x - μreg_Pi2)),

    # R:J_fATP
    J_fATP ~ E(μa_fATPt + μ_ATP_fe) - E(μa_fATPt + μ_ATP_fi),
    # R:J_mATP
    J_mATP ~ E(μa_mATPt + μ_ATP_me) - E(μa_mATPt + μ_ATP_mi),
    # R:J_fADP
    J_fADP ~ E(μa_fADPt + μ_ADP_fe) - E(μa_fADPt + μ_ADP_fi),
    # R:J_mADP
    J_mADP ~ E(μa_mADPt + μ_ADP_me) - E(μa_mADPt + μ_ADP_mi),
    # R:J_AMP
    J_AMP ~ E(μa_AMPt + μ_AMP_e) - E(μa_AMPt + μ_AMP_i),
    # R:J_Pi2
    J_Pi2 ~ E(μa_Pit + μ_Pi_e) - E(μa_Pit + μ_Pi_i),

    #R: AKi
    J_AKi ~ E(μa_AK + 2μ_ADP_fi) - E(μa_AK + μ_ATP_fi + μ_AMP_i),

    # R:ANT
    J_ANT ~ E(μa_ANT)*(1/(1+E(μ_ATP_fi-μ_ADP_fi+μ0_fADP-μ0_fATP-F*Psi_i)) - 1/(1+E(μ_ATP_fx-μ_ADP_fx+μ0_fADP-μ0_fATP-F*Psi_x)))/(1+E(μreg_ANT-μ_ADP_fi)),

    # R:J_Pi1
    J_Pi1 ~ (E(μa_PiHt + μ_Pi_i + 2μ_H_i) - E(μa_PiHt + μ_Pi_x + 2μ_H_x)) * 1/(1 + E(μ_H_i - μd_HPi + μ_Pi_i - μreg_PiHt)/(1 + E(μ_H_i - μd_HPi))) / (1 + E(μ_H_i - μd_HPi)),

    # C:Psi
    D(dPsi) ~ (2J_C3 + 4J_C1 + 4J_C4 - J_ANT - J_Hle - J_K - n_A*J_F1) / C_im,
    
    # Algebraic expressions
    dG_H ~ F*dPsi + μ_H_i - μ_H_x,
    Psi_x ~ -0.65dPsi,
    Psi_i ~ 0.35dPsi,
    W_x ~ matrix_mito_ratio()*W_m,
    W_i ~ ims_mito_ratio()*W_m,

    # Useful observables to plot
    ATP_x ~ ATP_fx + ATP_mx,
    ADP_x ~ ADP_fx + ADP_mx,
    ATP_i ~ ATP_fi + ATP_mi,
    ADP_i ~ ADP_fi + ADP_mi,
    ΔG_C1 ~ μ_NAD_x - (μ_NADH_x + μ_Q - μ_QH2 + μ_H_x - 4dG_H),
    ΔG_C3 ~ 2μ_Cred - (2μ_Cox - 4dG_H - 2μ_H_x + 2F*dPsi + μ_QH2 - μ_Q),
    ΔG_C4 ~ 2F*dPsi + 2μ_Cox + μ_H2O - (2μ_H_x + 2μ_Cred + 1/2*μ_O2 - 2dG_H),
    ΔG_F1 ~ μ_ATP_mx + μ_H2O - (μ_ADP_mx + μ_Pi_x + μ_H_x + n_A*dG_H),
    ΔG_Pi1 ~ μ_Pi_x + 2μ_H_x - (μ_Pi_i + 2μ_H_i)
]

"""
Differential equations and chemical potentials for additional species in the intermembrane space in the in vitro model.
"""
ims_species_invitro() = [
    # Ce:ATP_mi
    D(ATP_mi) ~ (J_MgATPi + J_mATP) / W_i,
    μ_ATP_mi ~ μ(μ0_mATP,ATP_mi),

    # Ce:ADP_mi
    D(ADP_mi) ~ (J_MgADPi + J_mADP) / W_i,
    μ_ADP_mi ~ μ(μ0_mADP,ADP_mi)
]


"""
Equations for the cytosolic species under experimental conditions (for comparison with the original Beard model).
"""
cytosolic_species_invitro() = Equation[
    # Se:Mg_i
    Mg_e ~ Mg_tot - ADP_me, # Constant concentration since ADP_me is constant
    Mg_i ~ Mg_e,
    μ_Mg_i ~ μ(μ0_Mg,Mg_i),
    
    # Se:ATP_me
    ATP_me ~ 0.5ATP_e + 0.5K_DT + 0.5Mg_tot - 0.5sqrt((ATP_e + K_DT + Mg_tot)^2.0 - 4.0ATP_e*Mg_tot), # Note: This assumes that ATP is set to zero. If ATP is nonzero and Mg_tot is used to calculate MgATP and MgADP, this equation needs to be updated.
    μ_ATP_me ~ μ(μ0_mATP,ATP_me),

    # Se:ATP_fe
    μ_ATP_fe ~ μ(μ0_fATP,ATP_fe),
    ATP_fe ~ ATP_e - ATP_me,

    # Se:ADP_me
    μ_ADP_me ~ μ(μ0_mADP,ADP_me),
    ADP_me ~ 0.5ADP_e + 0.5K_DD + 0.5Mg_tot - 0.5sqrt((ADP_e + K_DD + Mg_tot)^2.0 - 4.0ADP_e*Mg_tot),

    # Se:ADP_fe
    μ_ADP_fe ~ μ(μ0_fADP,ADP_fe),
    ADP_fe ~ ADP_e - ADP_me,

    # Se:μ_AMP_e
    μ_AMP_e ~ μ(μ0_AMP,AMP_e),

    # Se:Pi_e
    μ_Pi_e ~ μ(μ0_Pi,Pi_e)
]

"""
Equations for the bond graph model under in vitro conditions
"""
equations_bg_invitro() = [
    equations_bg_base();
    cytosolic_species_invitro();
    ims_species_invitro()
]

"""
Parameters for the bond graph model under in vitro conditions
"""
parameters_bg_invitro() = merge(
    parameters_beard(),
    bg_species_parameters(),
    bg_reaction_parameters(),
    bg_regulation_parameters(),
    buffering_parameters()
)

"""
States for the bond graph model under in vitro conditions
"""
sts_bg_invitro() = [
    dPsi, H_x, Pi_x, QH2, Cred, 
    AMP_i, ADP_fi, ATP_fx, ATP_mx, 
    Mg_x, NAD_x, NADH_x, ATP_fi, ATP_mi, 
    K_x, O2, ADP_mi, ADP_mx, Pi_i, 
    ADP_fx, B_x, HB_x
]

"""
Initial conditions for the in vitro bond graph model. Small concentrations are used to avoid numerical issues with taking the log of zero.
"""
function u0_bg_invitro()
    u0 = u0_beard()
    pb = parameters_beard()
    ADP_x0 = u0[ADP_x]
    u0[NAD_x] = pb[NADtot] - u0[NADH_x]

    delete!(u0, ADP_i)
    u0[ADP_fi] = 1e-15
    delete!(u0, ATP_x)
    u0[ATP_fx] = 1e-15
    delete!(u0, ATP_i)
    u0[ATP_fi] = 1e-15
    delete!(u0, ADP_x)
    u0[ADP_fx] = ADP_x0 - u0[ADP_mx]

    p_buff = buffering_parameters()
    u0[B_x] = p_buff[Btot]/(1+u0[H_x]/p_buff[K_dHbuff])
    u0[HB_x] = p_buff[Btot] - u0[B_x]

    return u0
end

"""
Returns a system of equations for the in vitro bond graph model.
"""
function bg_model_invitro()
    eqns = equations_bg_invitro()
    parameter_vals = parameters_bg_invitro()
    u0 = u0_bg_invitro()

    ps = [
        Qtot, AMP_e, K_DT, r, n_A, C_im, NADtot, F, RT, K_DT, Ctot, K_DD, Mg_tot, pH_e, W_m, Pi_e, ATP_e, ADP_e, K_e, k_mADP,
        μ0_NADH_x, μ0_NAD_x, μ0_Q, μ0_QH2, μ0_H, μ0_Cox, μ0_Cred, μ0_O2, μ_H2O, μ0_fATP, μ0_fADP, μ0_Pi, μ0_Mg, μ0_mATP, μ0_mADP, μ0_K, μ0_AMP, μ0_B, μ0_HB,
        μa_C1, μa_C3, μa_C4, μa_F1, μa_MgATP, μa_MgADP, μa_Hle, μa_KH, μa_DH, μ_r, μa_fATPt, μa_mATPt, μa_fADPt, μa_mADPt, μa_AMPt, μa_Pit, μa_AK, μa_ANT, μa_PiHt, μa_Hb,
        μreg_Pi1, μreg_Pi2, μreg_Pi3, μreg_Pi4, μreg_O2, μreg_Cred, μd_HPi, μreg_PiHt, μreg_ANT,
        Btot, K_dHbuff
    ]
    sts = sts_bg_invitro()
    
    odesys = ODESystem(eqns,t,sts,ps;name=:oxphos,defaults=merge(parameter_vals,u0))
    sys = structural_simplify(odesys)
    return sys
end

# Variables associated with creatine kinase
@parameters μa_ATPase μa_CK μ0_Cr μ0_PCr
@variables Mg_e_d(t) Mg_i_d(t) ATP_me_d(t) ATP_fe_d(t) ADP_me_d(t) ADP_fe_d(t) AMP_e_d(t) Pi_e_d(t)
@variables W_e(t) J_MgATPe(t) J_MgADPe(t) J_ATPase(t) J_AKe(t) J_CK(t) μ_Cr(t) μ_PCr(t) Cr(t) PCr(t)


"""
Returns the equations for cytosolic Mg, ATP and ADP.
The equations include additional reactions for buffering and ATP consumption.
A creatine kinase reaction has been added, as well as the transport of creatine.
The species option allows the volume fraction to be adjusted for rats or humans.
"""
function cytosolic_species_invivo(;species=:human,ATPase_rate=:constant)
    rates = Equation[
        D(Mg_e_d) ~ (-J_MgATPe-J_MgADPe-J_MgATPi-J_MgADPi)/(W_e+W_i),
        D(ATP_me_d) ~ (-J_mATP+J_MgATPe-J_ATPase-J_CK)/W_e,
        D(ATP_fe_d) ~ (-J_fATP-J_MgATPe+J_AKe)/W_e,
        D(ADP_me_d) ~ (-J_mADP+J_MgADPe+J_ATPase+J_CK)/W_e,
        D(ADP_fe_d) ~ (-J_fADP-J_MgADPe-2J_AKe)/W_e,
        D(AMP_e_d) ~ (-J_AMP+J_AKe)/W_e,
        D(Pi_e_d) ~ (-J_Pi2+J_ATPase)/W_e,
        D(Cr) ~ (-J_CK-J_Cr)/W_e,
        D(PCr) ~ (J_CK-J_PCr)/W_e,
        W_e ~ cyto_mito_ratio(species)*W_m 
    ]

    potentials = Equation[    
        # Ce:Mg_e
        Mg_i_d ~ Mg_e_d,
        μ_Mg_i ~ μ(μ0_Mg,Mg_i_d),
        # Ce:ATP_me
        μ_ATP_me ~ μ(μ0_mATP,ATP_me_d),
        # Ce:ATP_fe
        μ_ATP_fe ~ μ(μ0_fATP,ATP_fe_d),
        # Ce:ADP_me
        μ_ADP_me ~ μ(μ0_mADP,ADP_me_d),
        # Ce:ADP_fe
        μ_ADP_fe ~ μ(μ0_fADP,ADP_fe_d),
        # Ce:AMP_d
        μ_AMP_e ~ μ(μ0_AMP,AMP_e_d),
        # Ce:Pi_e
        μ_Pi_e ~ μ(μ0_Pi,Pi_e_d),
        # Ce:Cr
        μ_Cr ~ μ(μ0_Cr,Cr),
        # Ce:PCr
        μ_PCr ~ μ(μ0_PCr,PCr)
    ]

    if ATPase_rate == :constant
        eq_ATPase = (J_ATPase ~ E(μa_ATPase + μ_ATP_me + μ_H2O) - E(μa_ATPase + μ_ADP_me + μ_Pi_e + μ_H_i))
    else
        eq_ATPase = (J_ATPase ~ ATPase_rate)
    end
    fluxes = [
        # R:J_MgATPe
        J_MgATPe ~ E(μa_MgATP + μ_Mg_i + μ_ATP_fe) - E(μa_MgATP + μ_ATP_me),
        # R:J_MgADPe
        J_MgADPe ~ E(μa_MgADP + μ_Mg_i + μ_ADP_fe) - E(μa_MgADP + μ_ADP_me),
        #R:ATPase
        eq_ATPase,
        # R: AKe
        J_AKe ~ E(μa_AK + 2μ_ADP_fe) - E(μa_AK + μ_ATP_fe + μ_AMP_e),
        # R: CK
        J_CK ~ E(μa_CK + μ_ATP_me + μ_Cr) - E(μa_CK + μ_ADP_me + μ_PCr + μ_H_i)
    ]

    return [rates;potentials;fluxes]
end

@parameters μa_Cr μa_PCr μa_CKi
@variables Cr_i(t) PCr_i(t) μ_Cr_i(t) μ_PCr_i(t) J_CK_i(t) J_Cr(t) J_PCr(t)
"""
Retuns the additional equations for the intermembrane compartment in the human model.
"""
function ims_species_invivo()
    species = [
        # Ce: ATP_mi
        D(ATP_mi) ~ (J_MgATPi + J_mATP - J_CK_i) / W_i,
        μ_ATP_mi ~ μ(μ0_mATP,ATP_mi),

        # Ce: ADP_mi
        D(ADP_mi) ~ (J_MgADPi + J_mADP + J_CK_i) / W_i,
        μ_ADP_mi ~ μ(μ0_mADP,ADP_mi),

        # Ce: Cr_i
        D(Cr_i) ~ (-J_CK_i + J_Cr) / W_i,
        μ_Cr_i ~ μ(μ0_Cr,Cr_i),

        # Ce: PCr_i
        D(PCr_i) ~ (J_CK_i + J_PCr) / W_i,
        μ_PCr_i ~ μ(μ0_PCr,PCr_i)
    ]

    fluxes = [
        # R: CK_i
        J_CK_i ~ E(μa_CKi + μ_ATP_mi + μ_Cr_i) - E(μa_CKi + μ_ADP_mi + μ_PCr_i + μ_H_i),
        # R: J_Cr
        J_Cr ~ E(μa_Cr + μ_Cr) - E(μa_Cr + μ_Cr_i),
        # R: J_PCr
        J_PCr ~ E(μa_PCr + μ_PCr) - E(μa_PCr + μ_PCr_i)
    ]

    return [species;fluxes]
end

"""
Returns a dict with the initial conditions for the human model.

Notes:
Total ATP, total ADP, Pi, Cr and PCr concentrations in the myofibril were taken from the first row of Table 2 in Aliev and Saks (1997; https://doi.org/10.1016/S0006-3495(97)78082-2). The free and Mg-bound concentrations of ATP and ADP were calculated assuming equilibrium of Mg binding.

A free Mg concentration of 0.91mM was assumed, from Watanabe and Konishi, (2001; https://doi.org/10.1007/s004240000499). Initial concentrations in all compartments were initialised to this value.

ATP and ADP concentrations in the matrix were calculated using equation 2 of Aliev and Saks.

The AMP concentration was calculated assuming equilibrium with free ATP and ADP via the adenylate kinase reaction 2ADP ⇌ ATP + AMP.

Buffer concentrations initialised such that buffer reaction is at equilibrium, with total buffer concentration obtained from optimisation.
"""
function u0_bg_invivo()
    pb = parameters_beard()
    pbg = parameters_bg_invivo()
    u0_vitro = u0_bg_invitro()

    # Set concentrations from Aliev and Saks (https://doi.org/10.1016/S0006-3495(97)78082-2). Assume that the initial myofibril concentration is the same as the inner mitochondrial membrane.
    ATP0 = 9.0e-3
    ADP0 = 16.2e-6
    Pi0 = 4.8e-3
    Cr0 = 6.8e-3
    PCr0 = 16.2e-3

    # The free intracellular Mg concentration has been estimated to be 0.91mM (Watanabe and Konishi, 2001; https://doi.org/10.1007/s004240000499)
    Mg_f0 = 0.91e-3

    # Calculate free and Mg-bound ATP assuming equilibrium of the reaction Mg + ATP ⇌ MgATP, and similar for ADP
    fATP0 = ATP0 / (1 + Mg_f0/pb[K_DT])
    mATP0 = ATP0 - fATP0
    fADP0 = ADP0 / (1 + Mg_f0/pb[K_DD])
    mADP0 = ADP0 - fADP0

    # Calculate AMP concentration assuming equilibrium of the ADK reaction
    AMP0 = (fADP0^2/fATP0) * E(2pbg[μ0_fADP]-pbg[μ0_fATP]-pbg[μ0_AMP])

    # Use Equation 2 from Aliev and Saks to calculate matrix concentrations
    Mg_fx0 = Mg_f0
    ANPx0 = 10e-3
    ADP_x0 = ADP0 * ANPx0 / (ADP0 + ATP0/25)
    ATP_x0 = ANPx0 - ADP_x0

    (ATP_fx0,ATP_mx0) = free_and_bound_states(ATP_x0,Mg_fx0,pb[K_DT])
    (ADP_fx0,ADP_mx0) = free_and_bound_states(ADP_x0,Mg_fx0,pb[K_DD])

    Pi_x0 = 3e-3

    u0_cytosolic = Dict(
        Mg_e_d => Mg_f0,
        ATP_me_d => mATP0,
        ATP_fe_d => fATP0,
        ADP_me_d => mADP0,
        ADP_fe_d => fADP0,
        ATP_mi => mATP0,
        ATP_fi => fATP0,
        ADP_mi => mADP0,
        ADP_fi => fADP0,
        AMP_e_d => AMP0,
        AMP_i => AMP0,
        Pi_e_d => Pi0,
        Cr => Cr0,
        PCr => PCr0,
        Cr_i => Cr0,
        PCr_i => PCr0,
        Mg_x => Mg_fx0,
        ATP_fx => ATP_fx0,
        ATP_mx => ATP_mx0,
        ADP_fx => ADP_fx0,
        ADP_mx => ADP_mx0,
        Pi_x => Pi_x0
    )

    return merge(u0_vitro,u0_cytosolic)
end 

function free_and_bound_states(xtot,Mg,Kd)
    free = xtot/(1+Mg/Kd) 
    bound = xtot-free
    return (free,bound)
end

"""
States for the in vivo model
"""
sts_invivo() = [sts_bg_invitro(); Mg_e_d; ATP_me_d; ATP_fe_d; ADP_me_d; ADP_fe_d; AMP_e_d; Pi_e_d; Cr; PCr; Cr_i; PCr_i]

"""
Returns a dict of parameters for ATP consumption. Value calculated assuming mass action to achieve a reaction rate equivalent to ATP generation at standard conditions for the Beard model.
"""
parameters_ATPase() = Dict(μa_ATPase => 3487.09)

"""
Returns a dict of parameters for creatine kinase. The value of the reaction rate calculated assuming the reaction operates near equilibrium. The parameter for Cr was taken from Wu et al. (2007; https://doi.org/10.1074/jbc.M701024200), and the parameter for PCr was calculated assuming that the reaction PCr + ADP + H ⇌ Cr + ATP has an equilibrium constant of 3.57e8 M^-1 (as per Wu et al.). As a sanity check, a similar value (4.05e8 M^-1) was derived from Table 4 of Teague and Dobson (1992; https://doi.org/10.1016/S0021-9258(19)49682-8).

Currently, we assume that creatine is transported at approximately similar rates to inorganic phoshpate.
"""
function parameters_CK(species=:human)
    pbg = parameters_bg_invitro()
    pb = parameters_beard()
    μ0_Cr_val = -252.68
    μ0_PCr_val = μ0_Cr_val+pbg[μ0_mATP]-pbg[μ0_mADP]-pbg[μ0_H]+L(3.57e8*pb[K_DD]/pb[K_DT])
    μa_CK_val = L(1e6)-μ0_Cr_val-pbg[μ0_mATP]
    Dict(
        μa_CK => μa_CK_val,
        μa_CKi => μa_CK_val + L(ims_mito_ratio()/cyto_mito_ratio(species)),
        μa_Cr => L(300)-μ0_Cr_val,
        μa_PCr => L(300)-μ0_PCr_val,
        μ0_Cr => μ0_Cr_val,
        μ0_PCr => μ0_PCr_val
    )
end

"""
Returns parameters for the in vivo model
"""
parameters_bg_invivo(species=:human) = merge(parameters_bg_invitro(),parameters_ATPase(),parameters_CK(species))

"""
Returns a system of equations for the in vivo bond graph model.
"""
function bg_model_invivo(;stat=[],transformations=Dict(),ATPase_rate=:constant,species=:human)
    eqns = [
        equations_bg_base(); 
        cytosolic_species_invivo(ATPase_rate=ATPase_rate,species=species); 
        ims_species_invivo()
    ]
    for x in stat
        stat_var!(eqns,x)
    end
    eqns_new = substitute(eqns,transformations)
    p_vals = parameters_bg_invivo(species)
    odesys = ODESystem(eqns_new;name=:oxphos,defaults=merge(p_vals,u0_bg_invivo()))
    sys = structural_simplify(odesys)
    return sys
end

"""
Modifies equations so that the variable given by x is fixed over the duration of a simulation
"""
function stat_var!(eqns,x)
    for (i,eqn) in enumerate(eqns)
        if isequal(eqn.lhs,D(x))
            eqns[i] = eqn.lhs ~ 0
            break
        end
    end
end

"""
Load in data for PCr/ATP ratio
"""
function load_PCrATP_data()
    json = JSON.parsefile("../data/functional/Vendelin2000_PCrATP.json")
    dataset = json["datasetColl"]

    data_PCrATP = dataset[1]["data"]
    VO2_PCrATP = [p["value"][1] for p in data_PCrATP]
    PCrATP_data = [p["value"][2] for p in data_PCrATP]

    data_PCrCr = dataset[2]["data"]
    VO2_PCrCr = [p["value"][1] for p in data_PCrCr]
    PCrCr_data = [p["value"][2] for p in data_PCrCr]
    
    return (VO2_PCrATP,PCrATP_data)
end

"""
Returns an array of activity parameters for ATPase, for use in fitting.
"""
function ATPase_fit_params()
    pbg = parameters_bg_invivo()
    ATP_consumption_ratio = LinRange(0.01,7,20)
    ATPase_params = @. pbg[μa_ATPase] + L(ATP_consumption_ratio)
end

"""
Simulate the oxidative phosphorylation model for fitting purposes. For each workload (determined by ATPase activity), the rate of complex IV and PCr/ATP ratio are collected. The flux through complex IV is converted to VO2 (in μmol/min/g dry wt) using the conversion factor in Ghosh et al. (2018; https://doi.org/10.1371/journal.pcbi.1006640).
"""
function simulate_VO2(sys,μa_ANT_val;ATPase_params=ATPase_fit_params(),tspan=(0.0,100.0),alg=CVODE_BDF())
    ps = [[μa_ATPase => m, μa_ANT => μa_ANT_val] for m in ATPase_params]
    sols = multi_sim(sys,tspan,ps;alg=alg)
    O2conv_vendelin = 54921 # From Ghosh et al. (2018; https://doi.org/10.1371/journal.pcbi.1006640)
    VO2 = O2conv_vendelin*extract_var(sols,J_C4)
    PCr_sol = extract_var(sols,PCr)
    ATP_sol = extract_var(sols,ATP_fe_d+ATP_me_d)
    PCrATP = PCr_sol./ATP_sol
    return (VO2, PCrATP)
end

"""
Returns a cost function for fitting the μa_ANT parameter. The cost function includes an error term to account for deviation from data and a regularisation term to account for deviations from a prior value determined from the Beard model.
"""
function PCr_cost(sys,μa_ANT_val; data=load_PCrATP_data())
    (VO2_PCrATP,PCrATP_data) = data
    pbg = parameters_bg_invivo(:rat)
    (x,y) = simulate_VO2(sys,μa_ANT_val)
    interp = LinearInterpolation(x,y,extrapolation_bc = NaN)
    predictions = interp.(VO2_PCrATP)
    return sum(abs2, PCrATP_data-predictions) + 0.01*(μa_ANT_val-pbg[μa_ANT])^2
end

"""
Fit the parameter for ANT using an exhaustive search.
"""
function fit_ANT_to_PCrATP()
    pbg = parameters_bg_invivo(:rat)
    μa_ANT_vals = pbg[μa_ANT] .+ LinRange(L(0.5),L(100),100)
    sys = bg_model_invivo(species=:rat)
    tspan = (0.0,100.0)
    options = Dict(:reltol => 1e-7, :abstol => 1e-10, :alg => CVODE_BDF())
    data = data=load_PCrATP_data()
    costs = [PCr_cost(sys,m;data=data) for m in μa_ANT_vals]
    copt = minimum(costs[costs .!== NaN])
    popt = μa_ANT_vals[findfirst(costs .== copt)]
    return popt,(μa_ANT_vals,costs)
end

"""
Returns the fitted value of μa_ANT. run_fit can be set to true to run the optimisation problem again.
"""
updated_μa_ANT(;run_fit=false) = run_fit ? fit_ANT_to_PCrATP() : -8.121878217633974

end