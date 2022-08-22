include("interfaz.jl")

function timeevolution(u, v, p, ρ, E, Δξ, Δη, γ, ξ, η, id)
    ρu_px, ρu_mx, ρu_py, ρu_my = minmood(ρ .*u, Δξ, Δη)
    ρv_px, ρv_mx, ρv_py, ρv_my = minmood(ρ .*v, Δξ, Δη)
    ρ_px, ρ_mx, ρ_py, ρ_my = minmood(ρ, Δξ, Δη)
    E_px, E_mx, E_py, E_my = minmood(E, Δξ, Δη)
    # -----------------------------------------------------------------------
    u_px, u_mx, u_py, u_my = ρu_px./ρ_px, ρu_mx./ρ_mx, ρu_py./ρ_py, ρu_my./ρ_my
    v_px, v_mx, v_py, v_my = ρv_px./ρ_px, ρv_mx./ρ_mx, ρv_py./ρ_py, ρv_my./ρ_my
    p_px = (E_px .- 0.5*ρ_px .*(u_px.^2 .+ v_px.^2))*(γ-1)
    p_mx = (E_mx .- 0.5*ρ_mx .*(u_mx.^2 .+ v_mx.^2))*(γ-1)
    p_py = (E_py .- 0.5 *ρ_py .*(u_py.^2 .+ v_py.^2))*(γ-1)
    p_my = (E_my .- 0.5 *ρ_my .*(u_my.^2 .+ v_my.^2))*(γ-1)
    # ----------------------------------------------------------------------
    s_px = soundspeed(γ, p_px, ρ_px)
    s_mx = soundspeed(γ, p_mx, ρ_mx)
    s_py = soundspeed(γ, p_py, ρ_py)
    s_my = soundspeed(γ, p_my, ρ_my)
    a_max = maxlocspeeds(u_px, u_mx, s_px, s_mx)[3]
    b_max = maxlocspeeds(v_py, v_my, s_py, s_my)[3]
    n = size(u, 1); m = size(u, 2)
    ξ_int, η_int = interfaz(ξ, η)
    h_ηint = metrics(ξ_int, η, id)[2]
    h_ξint = metrics(ξ, η_int, id)[1]
    a_h = h_ηint .*a_max
    b_h =  h_ξint .*b_max
    h_ξ, h_η, dhξ_dη, dhη_dξ = metrics(ξ, η, id)
    ΔV = (Δξ*Δη).*h_ξ .*h_η
    acoc = zeros(n, m)
    bcoc = zeros(n, m)
    for i=2:n-1
        for j=2:m-1
        acoc[i, j] = Δη * max(a_h[i, j], a_h[i-1,j]) /ΔV[i, j]
        bcoc[i, j] = Δξ  * max(b_h[i, j],b_h[i, j-1]) /ΔV[i, j]
        end
    end
    #a_Δξ = zeros(length(a_max))
    #b_Δη = zeros(length(b_max))
    Δt = (1/8)/(maximum([acoc bcoc]))
    return Δt
end

#function RHS(u, v, p, ρ, E, ξ, η, Δξ, Δη, id, signFp, signFm,  signGp, signGm)
function RHS(com, flux_um, flux_up, vel_p, vel_m, U_p, U_m, signR, signL, singP)
    flux = com.*(signR*vel_p .*flux_um .- signL*vel_m .*flux_up .-
                 singP*vel_p .*vel_m .*(U_m .- U_p))
    #G = com2.*(signGm*b_p .*G_um .- signGp*b_m .*G_up .-
    #            b_p .*b_m .*(ρ_my .- ρ_py))
    #F2 = com1.*(signFm*a_p .*F2_um .- signFp*a_m .*F2_up .-
    #            a_p .*a_m .*(ρ_mx .*u_mx .- ρ_px .*u_px))
    #G2 = com2.*(signGm*b_p .*G2_um .- signGp*b_m .*G2_up .-
    #            b_p .*b_m .*(ρ_my .*u_my .- ρ_py .*u_py))
    #F3 = com1.*(signFm*a_p .*F3_um .- signFp*a_m .*F3_up .-
    #            a_p .*a_m .*(ρ_mx .*v_mx .- ρ_px .*v_px))
    #G3 = com2.*(signGm*b_p .*G3_um .- signGp*b_m .*G3_up .-
    #            b_p .*b_m .*(ρ_my .*v_my .- ρ_py .*v_py))
    #F4 = com1.*(signFm*a_p .*F1_um .- signFp*a_m .*F1_up .-
    #            a_p .*a_m .*(E_mx .- E_px))
    #G4 = com2.*(signGm*b_p .*G1_um .- signGp*b_m .*G1_up .-
    #             b_p .*b_m .*(E_my .- E_py))
    #ΔV = (Δξ*Δη).*h_ξ .*h_η
    return flux
end




function func!(dU, U, par, t)
    γ = par[1]
    ξ = par[2]; η = par[3]
    Δξ = par[4]; Δη = par[5]
    id = par[6]; B_C = par[7]
    # -------------------------------------------------------------------
    ρ_px, ρ_mx, ρ_py, ρ_my = minmood(U.ρ, Δξ, Δη)
    ρu_px, ρu_mx, ρu_py, ρu_my = minmood(U.ρu, Δξ, Δη)
    ρv_px, ρv_mx, ρv_py, ρv_my = minmood(U.ρv, Δξ, Δη)
    E_px, E_mx, E_py, E_my = minmood(U.E, Δξ, Δη)
    # ------------------------------------------------------------------
    u_px, u_mx, u_py, u_my = ρu_px./ρ_px, ρu_mx./ρ_mx, ρu_py./ρ_py, ρu_my./ρ_my
    v_px, v_mx, v_py, v_my = ρv_px./ρ_px, ρv_mx./ρ_mx, ρv_py./ρ_py, ρv_my./ρ_my
    p_px = (E_px .- 0.5*ρ_px .*(u_px.^2 .+ v_px.^2))*(γ-1)
    p_mx = (E_mx .- 0.5*ρ_mx .*(u_mx.^2 .+ v_mx.^2))*(γ-1)
    p_py = (E_py .- 0.5*ρ_py .*(u_py.^2 .+ v_py.^2))*(γ-1)
    p_my = (E_my .- 0.5*ρ_my .*(u_my.^2 .+ v_my.^2))*(γ-1)
    # -------------------------------------------------------------------
    s_px = soundspeed(γ, p_px, ρ_px)
    s_mx = soundspeed(γ, p_mx, ρ_mx)
    s_py = soundspeed(γ, p_py, ρ_py)
    s_my = soundspeed(γ, p_my, ρ_my)
    a_p, a_m, a_max = maxlocspeeds(u_px, u_mx, s_px, s_mx)
    b_p, b_m, b_max = maxlocspeeds(v_py, v_my, s_py, s_my)
    ξ_int, η_int = interfaz(ξ, η)
    h_ηint = metrics(ξ_int, η, id)[2]
    h_ξint = metrics(ξ, η_int, id)[1]
    h_ξ, h_η, dhξ_dη, dhη_dξ = metrics(ξ, η, id)
    # ----------------------------------------------------------------
    F1_up, F2_up, F3_up, F4_up = FluxF(u_px, v_px, ρ_px, p_px, E_px)
    F1_um, F2_um, F3_um, F4_um = FluxF(u_mx, v_mx, ρ_mx, p_mx, E_mx)
    G1_up, G2_up, G3_up, G4_up = FluxG(u_py, v_py, ρ_py, p_py, E_py)
    G1_um, G2_um, G3_um, G4_um = FluxG(u_my, v_my, ρ_my, p_my, E_my)
    # ----------------------------------------------------------------
    ρ = U.ρ; u = U.ρu ./ρ; v =  U.ρv ./ρ;
    p = (U.E .- 0.5*ρ .*(u.^2 .+ v.^2))*(γ-1)
    S1, S2 = source(u, v, ρ, p, h_ξ, h_η, dhξ_dη, dhη_dξ)
    com1 = (Δη .*h_ηint)./(a_p .+ a_p)
    com2 = (Δξ .*h_ξint)./(b_p .+ b_p)
    ΔV = (Δξ*Δη).*h_ξ .*h_η
    # ------------------------------------------------------------------------
    F1 = RHS(com1, F1_up, F1_um, a_p, -a_p, ρ_px, ρ_mx, 1, 1, 1)
    G1 = RHS(com2, G1_up, G1_um, b_p, -b_p, ρ_py, ρ_my, 1, 1, 1)
    F2 = RHS(com1, F2_up, F2_um, a_p, -a_p, ρ_px .*u_px, ρ_mx .*u_mx, 1, 1, 1)
    G2 = RHS(com2, G2_up, G2_um, b_p, -b_p, ρ_py .*u_py, ρ_my .*u_my, 1, 1, 1)
    F3 = RHS(com1, F3_up, F3_um, a_p, -a_p, ρ_px .*v_px, ρ_mx .*v_mx, 1, 1, 1)
    G3 = RHS(com2, G3_up, G3_um, b_p, -b_p, ρ_py .*v_py, ρ_my .*v_my, 1, 1, 1)
    F4 = RHS(com1, F4_up, F4_um, a_p, -a_p, E_px, E_mx, 1, 1, 1)
    G4 = RHS(com2, G4_up, G4_um, b_p, -b_p, E_py, E_my, 1, 1, 1)

    # ---------------------------------------------------------
    n = size(u, 1); m = size(u, 2)
    for i=2:n-1
        for j=2:m-1
            dU.ρ[i, j] = (1/ΔV[i,j])*(F1[i-1,j] - F1[i,j] + G1[i,j-1] - G1[i,j])
            dU.ρu[i, j] = (1/ΔV[i,j])*(F2[i-1,j] - F2[i,j] + G2[i,j-1] - G2[i,j]) + S1[i,j]
            dU.ρv[i, j] = (1/ΔV[i,j])*(F3[i-1,j] - F3[i,j] + G3[i,j-1] - G3[i,j]) + S2[i,j]
            dU.E[i, j] = (1/ΔV[i,j])*(F4[i-1,j] - F4[i,j] + G4[i,j-1] - G4[i,j])
        end
    end
   for i=1:4
       if B_C[i, 2] == "Interfaz"
           pos = B_C[i, 1]
           if pos == 1
                # Determinació  coordenadas y étricas cuarto de celda
                ξ_quart = quarterint(ξ, 1, "x", Δξ)
                h_ξquar, h_ηquart, dhξ_dηquart, dhη_dξquart = metrics(ξ_quart, η, id)
                h_ξquarint = metrics(ξ_quart, η_int, id)[1]
                com2qua = ((Δξ/2) .*h_ξquarint)./(b_p .+ b_p)
                ΔV_1 = (Δξ/2)*Δη .*h_ξquar .*h_ηquart
                # Cálculo de flujos numéricos
                F11 = RHS(com1, F1_up, F1_um, a_p, -a_p, ρ_px, ρ_mx, 1, -1, -1)[1, :]
                F22 = RHS(com1, F2_up, F2_um, a_p, -a_p, ρ_px .*u_px, ρ_mx .*u_mx,
                         1, -1, -1)[1,:]
                F33 = RHS(com1, F3_up, F3_um, a_p, -a_p, ρ_px .*v_px, ρ_mx .*v_mx,
                         1, -1, -1)[1,:]
                F44 = RHS(com1, F4_up, F4_um, a_p, -a_p, E_px, E_mx, 1, -1, -1)[1,:]
                G11 = RHS(com2qua, G1_up, G1_um, b_p, -b_p, ρ_py, ρ_my, 1, 1, 1)[1,:]
                G22 = RHS(com2qua, G2_up, G2_um, b_p, -b_p, ρ_py .*u_py,
                          ρ_my .*u_my, 1, 1, 1)[1,:]
                G33 = RHS(com2qua, G3_up, G3_um, b_p, -b_p, ρ_py .*v_py,
                         ρ_my .*v_my, 1, 1, 1)[1,:]
                G44 = RHS(com2qua, G4_up, G4_um, b_p, -b_p, E_py, E_my, 1, 1, 1)[1,:]
                # Cálculo fuente
                S11, S22 = source(u, v, ρ, p, h_ξquar, h_ηquart,
                                  dhξ_dηquart, dhη_dξquart)
             for j=2:m-1
                dU.ρ[1, j] = (1/ΔV_1[1, j])*(F11[j] + (G11[j-1] - G11[j]))
                dU.ρu[1, j] = (1/ΔV_1[1, j])*(F22[j] + (G22[j-1] - G22[j])) + S11[1,j]
                dU.ρv[1, j] = (1/ΔV_1[1, j])*(F33[j] + (G33[j-1] - G33[j])) + S22[1,j]
                dU.E[1, j] = (1/ΔV_1[1, j])*(F44[j] + (G44[j-1] - G44[j]))
             end
          end
           if pos == 2
                 # Determinació  coordenadas y étricas cuarto de celda
                 ξ_quart = quarterint(ξ, -1, "x", Δξ)
                 h_ξquar, h_ηquart, dhξ_dηquart, dhη_dξquart = metrics(ξ_quart, η, id)
                 h_ξquarint = metrics(ξ_quart, η_int, id)[1]
                 com2qua = ((Δξ/2) .*h_ξquarint)./(b_p .+ b_p)
                 ΔV_1 = (Δξ/2)*Δη .*h_ξquar .*h_ηquart
                 F11 = RHS(com1, F1_up, F1_um, a_p, -a_p, ρ_px, ρ_mx, -1, 1, 1)[end,:]
                 F22 = RHS(com1, F2_up, F2_um, a_p, -a_p, ρ_px .*u_px, ρ_mx .*u_mx,
                         -1, 1, 1)[end,:]
                 F33 = RHS(com1, F3_up, F3_um, a_p, -a_p, ρ_px .*v_px, ρ_mx .*v_mx,
                         -1, 1, 1)[end,:]
                 F44 = RHS(com1, F4_up, F4_um, a_p, -a_p, E_px, E_mx, -1, 1, 1)[end,:]
                 G11 = RHS(com2qua, G1_up, G1_um, b_p, -b_p, ρ_py, ρ_my, 1, 1, 1)[end,:]
                 G22 = RHS(com2qua, G2_up, G2_um, b_p, -b_p, ρ_py .*u_py,
                           ρ_my .*u_my, 1, 1, 1)[end,:]
                 G33 = RHS(com2qua, G3_up, G3_um, b_p, -b_p, ρ_py .*v_py,
                          ρ_my .*v_my, 1, 1, 1)[end,:]
                 G44 = RHS(com2qua, G4_up, G4_um, b_p, -b_p, E_py, E_my, 1, 1, 1)[end,:]
                 # Cálculo fuente
                 S11, S22 = source(u, v, ρ, p, h_ξquar, h_ηquart,
                                   dhξ_dηquart, dhη_dξquart)
              for j=2:m-1
                 dU.ρ[end, j] = (1/ΔV_1[end, j])*(F11[j] + (G1[j-1] - G1[j]))
                 dU.ρu[end, j] = (1/ΔV_1[end, j])*(F22[j] + (G2[j-1] - G2[j])) + S1[end,j]
                 dU.ρv[end, j] = (1/ΔV_1[end, j])*(F33[j] + (G3[j-1] - G3[j])) + S2[end,j]
                 dU.E[end, j] = (1/ΔV_1[end, j])*(F44[j] + (G4[j-1] - G4[j]))
              end
            end
            if pos == 3
                   ind = 1
                   # Determinació  coordenadas y étricas cuarto de celda
                   η_quart = quarterint(η, 1, "y", Δη)
                   h_ξquar, h_ηquart, dhξ_dηquart, dhη_dξquart = metrics(ξ, η_quart, id)
                   h_ηquarint = metrics(ξ_int, η_quart, id)[2]
                   com1quar = ((Δη/2) .* h_ηquarint)./(a_p .+ a_p)
                   ΔV_1 = (Δη/2)*Δξ .*h_ξquar .*h_ηquart
                   # Determinaciín flujos numéricos
                   G1_1 = RHS(com2, G1_up, G1_um, b_p, -b_p, ρ_py, ρ_my, 1, -1, -1)[:,ind]
                   G2_1 = RHS(com2, G2_up, G2_um, b_p, -b_p, ρ_py .*u_py,
                            ρ_my .*u_my, 1, -1, -1)[:,ind]
                   G3_1 = RHS(com2, G3_up, G3_um, b_p, -b_p, ρ_py .*v_py,
                            ρ_my .*v_my, 1, -1, -1)[:,ind]
                   G4_1 = RHS(com2, G4_up, G4_um, b_p, -b_p, E_py, E_my, 1, -1, -1)[:,ind]
                   F1_1 = RHS(com1quar, F1_up, F1_um, a_p, -a_p, ρ_px, ρ_mx, 1, 1, 1)[:,ind]
                   F2_2 = RHS(com1quar, F2_up, F2_um, a_p, -a_p, ρ_px .*u_px,
                              ρ_mx .*u_mx, 1, 1, 1)[:,ind]
                   F3_3 = RHS(com1quar, F3_up, F3_um, a_p, -a_p, ρ_px .*v_px,
                              ρ_mx .*v_mx, 1, 1, 1)[:,ind]
                   F4_4 = RHS(com1quar, F4_up, F4_um, a_p, -a_p, E_px, E_mx, 1, 1, 1)[:,ind]
                   S11, S22 = source(u, v, ρ, p, h_ξquar, h_ηquart,
                                     dhξ_dηquart, dhη_dξquart)
                for i=2:n-1
                   dU.ρ[i, ind] = (1/ΔV_1[i, 1])*(G1_1[i] + (F1_1[i-1] - F1_1[i]))
                   dU.ρu[i, ind] = (1/ΔV_1[i, 1])*(G2_1[i] + (F2_2[i-1] - F2_2[i])) + S11[i,ind]
                   dU.ρv[i, ind] = (1/ΔV_1[i, 1])*(G3_1[i] + (F3_3[i-1] - F3_3[i])) + S22[i,ind]
                   dU.E[i, ind] = (1/ΔV_1[i, 1])*(G4_1[i] + (F4_4[i-1] - F4_4[i]))
                end
            end
            if pos == 4
                ind = m
                # Determinació  coordenadas y étricas cuarto de celda
                η_quart = quarterint(η, -1, "y", Δη)
                h_ξquar, h_ηquart, dhξ_dηquart, dhη_dξquart = metrics(ξ, η_quart, id)
                h_ηquarint = metrics(ξ_int, η_quart, id)[2]
                com1quar = ((Δη/2) .* h_ηquarint)./(a_p .+ a_p)
                ΔV_1 = (Δη/2)*Δξ .*h_ξquar .*h_ηquart

                # Determinaciín flujos numéricos
                G1_1 = RHS(com2, G1_up, G1_um, b_p, -b_p, ρ_py, ρ_my, -1, 1, 1)[:,end]
                G2_1 = RHS(com2, G2_up, G2_um, b_p, -b_p, ρ_py .*u_py,
                         ρ_my .*u_my, -1, 1, 1)[:,end]
                G3_1 = RHS(com2, G3_up, G3_um, b_p, -b_p, ρ_py .*v_py,
                         ρ_my .*v_my, -1, 1, 1)[:,end]
                G4_1 = RHS(com2, G4_up, G4_um, b_p, -b_p, E_py, E_my, -1, 1, 1)[:,end]
                F1_1 = RHS(com1quar, F1_up, F1_um, a_p, -a_p, ρ_px, ρ_mx, 1, 1, 1)[:,end]
                F2_2 = RHS(com1quar, F2_up, F2_um, a_p, -a_p, ρ_px .*u_px,
                           ρ_mx .*u_mx, 1, 1, 1)[:,end]
                F3_3 = RHS(com1quar, F3_up, F3_um, a_p, -a_p, ρ_px .*v_px,
                           ρ_mx .*v_mx, 1, 1, 1)[:,end]
                F4_4 = RHS(com1quar, F4_up, F4_um, a_p, -a_p, E_px, E_mx, 1, 1, 1)[:,end]
                S11, S22 = source(u, v, ρ, p, h_ξquar, h_ηquart,
                                  dhξ_dηquart, dhη_dξquart)
                for i=2:n-1
                   dU.ρ[i, ind] = (1/ΔV_1[i, end])*(G1_1[i] + (F1_1[i-1] - F1_1[i]))
                   dU.ρu[i, ind] = (1/ΔV_1[i, end])*(G2_1[i] + (F2_2[i-1] - F2_2[i])) + S11[i,end]
                   dU.ρv[i, ind] = (1/ΔV_1[i, end])*(G3_1[i] + (F3_3[i-1] - F3_3[i])) + S22[i,end]
                   dU.E[i, ind] = (1/ΔV_1[i, end])*(G4_1[i] + (F4_4[i-1] - F4_4[i]))
                end
            end
       end
   end
   return dU
end
