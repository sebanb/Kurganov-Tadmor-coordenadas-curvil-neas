##Ininital Conditions
function Initialconditions()
u_xini = 1038
u_yini = 0
p_ini = 101325.0
T_ini = 298.0
return u_xini, u_yini, p_ini, T_ini
end

function time()
    t0 = 0.000564588021779721
    tf = 0.06
    γ = 1.4
    R  = 287.03
    name = "Wedge refine5"
    return t0, tf, γ, R, name
end


function BoundaryCondition(Cor)
    df = []
#-------------------------------------------------------------------------------
# first sub-domain

θ = Cor[2]
co = cos.(θ); sn = sin.(θ)
left= Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Neumman",
                         "Neumman",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                          0,
                           0,
                           0]
               )
right = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Dirichlet",
                         "Dirichlet",
                         "Dirichlet",
                         "Dirichlet"],
               :value => [ -1038 *co,
                           1038 *sn,
                           101325.0,
                           298]
               )

top = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Interfaz",
                         "Interfaz",
                         "Interfaz",
                         "Interfaz"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )

bot = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Neumman",
                         "Dirichlet",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )

push!(df, [left, right, bot, top])
#-------------------------------------------------------------------------------
# second sub-domain


left= Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Interfaz",
                         "Interfaz",
                         "Interfaz",
                         "Interfaz"],
               :value => [ 0,
                          0,
                           0,
                           0]
               )

right = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Interfaz",
                         "Interfaz",
                         "Interfaz",
                         "Interfaz"],
               :value => [ 0,
                          0,
                           0,
                           0]
               )

bot = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => ["Neumman",
                         "Dirichlet",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )

top = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => ["Neumman",
                         "Neumman",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )

push!(df, [left, right, bot, top])
#-------------------------------------------------------------------------------
# Third sub-domain

left= Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Neumman",
                         "Neumman",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )

right = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Neumman",
                         "Neumman",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                          0,
                           0,
                           0]
               )

top = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => ["Neumman",
                         "Neumman",
                         "Neumman",
                         "Neumman"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )

bot = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => ["Interfaz",
                         "Interfaz",
                         "Interfaz",
                         "Interfaz"],
               :value => [ 0,
                           0,
                           0,
                           0]
               )
push!(df, [left, right, bot, top])
#--------------------------------------------------------------------------

return df
end
