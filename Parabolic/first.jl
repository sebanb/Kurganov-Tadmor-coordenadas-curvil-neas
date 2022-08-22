using DifferentialEquations
using ComponentArrays
using JLD
using WriteVTK

include("interfaz.jl")
include("meshes4.jl")
include("solver.jl")
include("BaIC4.jl")
include("Coordinatesbb.jl")



Dom, Cor, Δ, bc = Mesh()
t0, tf, γ, R, name = time()

if t0 == 0
    t = t0
    u_xini, u_yini, p_ini, T_ini = Initialconditions()
    u = []; v = [];p = [];E = [];ρ = []
    for  i=1:Dom[end,1]
        ξ = Cor[2*i-1]; η = Cor[2*i]
        n = size(Cor[2*i-1], 1)
        m = size(Cor[2*i], 1)
        u1 = fill!(zeros(n, m), 0)
        v1 = fill!(zeros(n, m), 0)
        E1 = fill!(zeros(n, m), 0)
        p1 = fill!(zeros(n, m), p_ini)
        T1 = fill!(zeros(n, m), T_ini)
        ρ1 = p1 ./ (R * T1)
        for j = 1:n
          for k = 1:m
           u1[j, k], v1[j, k] = Car2Cur(u_xini, u_yini, ξ[j], η[k], Dom[i, 2])
           E1[j, k] = ρ1[j, k] * (0.5 * (u1[j, k]^2 + v1[j, k]^ 2) +
                                 p1[j, k]/(ρ1[j, k] * (γ - 1)))
          end
        end
        push!(u, u1); push!(v, v1); push!(p, p1); push!(ρ, ρ1); push!(E, E1)
    end
else
    u = load(name*string(t0)*".jld")["u_ξ"]
    v = load(name*string(t0)*".jld")["u_η"]
    p = load(name*string(t0)*".jld")[ "p"]
    #E = load(name*".jld")["E"]
    #ρ = load(name*".jld")["ρ"]
    t = load(name*string(t0)*".jld")["t"]
    T = load(name*string(t0)*".jld")["T"]
end

BondC = BoundaryCondition(Cor)

t = t0
# Start the time interation
while t < tf
    Δt_dom = []
    for i = 1:Dom[end, 1]
        Δξ = Δ[2*i-1]
        Δη = Δ[2*i]
        Δt_sub = timeevolution(u[i], v[i], p[i], ρ[i], E[i], Δξ, Δη, γ)
        push!(Δt_dom, Δt_sub)
    end
    Δt = minimum(Δt_dom)
    tspan = (t, t + Δt)
    utot = []
    vtot = []
    Etot = []
    ptot = []
    Ttot = []
    for i = 1:Dom[end, 1]
        ξ = Cor[2*i-1]
        η = Cor[2*i]
        Δξ = Δ[2*i-1]
        Δη = Δ[2*i]
        id = Dom[i, 2]
        par = [γ, ξ, η, Δξ, Δη, id, bc[i]]
        U = ComponentArray(
            ρ = ρ[i],
            ρu = ρ[i] .* u[i],
            ρv = ρ[i] .* v[i],
            E = E[i],
        )
        prob = ODEProblem(func!, U, tspan, par)
        sol =
            solve(prob, RK4(), dt = Δt, save_start = false, adaptive = false)
        ρdom = sol.u[1].ρ
        udom = sol.u[1].ρu ./ ρdom
        vdom = sol.u[1].ρv ./ ρdom
        Edom = sol.u[1].E
        # dU = zero(U)
        # dU = func!(dU, U, par, t)
        # U  +=  Δt*dU
        # ρdom = U.ρ
        # udom = U.ρu ./ρdom
        # vdom = U.ρv ./ρdom
        # Edom = U.E
        pdom = (Edom .- 0.5 .* ρdom .* (udom .^ 2 .+ vdom .^ 2)) .* (γ - 1)
        Tdom = pdom ./ (R * ρdom)
        push!(utot, udom)
        push!(vtot, vdom)
        push!(ptot, pdom)
        push!(Ttot, Tdom)
    end

    global u = deepcopy(utot)
    global v = deepcopy(vtot)
    global p = deepcopy(ptot)
    global T = deepcopy(Ttot)
    for i = 1:Dom[end, 1]
        for j = 1:4
            B_c = bc[i][j, 2]
            Domn = bc[i][j, 3]
            Boundary_n = bc[i][j, 4]
            if B_c == "Interfaz"
                uint, row, col =
                    Interfaz(i, Domn, j, Boundary_n, utot, Δ, Cor, Dom, vtot)
                u[i][row, col] = uint
                #---------------
                vint, row, col =
                    Interfaz(i, Domn, j, Boundary_n, vtot, Δ, Cor, Dom, utot)
                v[i][row, col] = vint
                #---------------
                pint, row, col =
                    Interfaz(i, Domn, j, Boundary_n, ptot, Δ, Cor, Dom, ptot)
                p[i][row, col] = pint
                #---------------
                Tint, row, col =
                    Interfaz(i, Domn, j, Boundary_n, Ttot, Δ, Cor, Dom, Ttot)
                T[i][row, col] = Tint
            end
        end
    end
    # Boundary Condition
    dictioanry_var = Dict("u_ξ" => u, "u_η" => v, "p" => p, "T" => T)
    #u = deepcopy(utot)
    for i = 1:Dom[end, 1]
        for j = 1:4#4#4
            for v = 1:4#4
                var = get(dictioanry_var, BondC[i][j][:var][v], 0)
                Cond = BondC[i][j][:Cond][v]
                Value = BondC[i][j][:value][v]
                var[i] = BC(i, j, var, Cond, Value)
            end
        end
        global ρ[i] = p[i] ./ (R * T[i])
        global E[i] =
            ρ[i] .* (0.5 * (u[i] .^ 2 .+ v[i] .^ 2) .+ p[i] ./ (ρ[i] * (γ - 1)))
    end

    global t = t + Δt
end

n = length(u[1])
save(name*"nelem="*string(n)*".jld", "t", t, "u_ξ", u, "u_η", v, "p", p, "T", T, "ρ", ρ, "E", E)
x = []; y = [];
for i = 1:Dom[end, 1]
    if i == 1
    ξ = - Cor[2*i-1]
    η = - Cor[2*i]
    else
    ξ = Cor[2*i-1]
     η = Cor[2*i]
    end

    n = size(Cor[2*i-1], 1)
    m = size(Cor[2*i], 1)
    x1 = zeros(n, m); y1 = zeros(n, m)
    u1 = zeros(n, m); v1 = zeros(n, m)
    for j=1:n
        for k=1:m
            x1[j, k], y1[j, k] = changcoor(ξ[j], η[k],Dom[i, 2])
        end
    end
    push!(x, x1); push!(y, y1);
end

save(name*"nelem="*string(n)*".jld", "t", t, "u_ξ", u, "u_η", v, "p", p, "T", T, "ρ", ρ, "E", E, "y", y)
saved_files = vtk_multiblock("full_domain_nelem="*string(n)) do vtm
 for i=1:Dom[end, 1]
    vtk_grid(vtm, x[i], y[i]) do vtk
        vtk["Pressure"] = p[i]
        vtk["Temperature"] = T[i]
        vtk["u_ξ"] = u[i]
        vtk["u_η"] = v[i]
     end
 end
end
