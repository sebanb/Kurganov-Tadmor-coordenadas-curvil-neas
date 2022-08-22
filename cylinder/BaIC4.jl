#Ininital Conditions
function Initialconditions()
u_xini = 692
u_yini = 0
p_ini = 101325.0
T_ini = 298.0
return u_xini, u_yini, p_ini, T_ini
end

function time()
    t0 = 0
    tf = 0.3
    γ = 1.4
    R  = 287.03
    name = "BluntBody"
    return t0, tf, γ, R, name
end


function BoundaryCondition(Cor)
    df = []
#-------------------------------------------------------------------------------
# first sub-domain
θ = Cor[2]
co = cos.(θ); sn = sin.(θ)



left= Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Dirichlet",
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
               :value => [ -692 .*co,
                          692 .*sn,
                          101325,
                          298]
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

top = Dict(:var => ["u_ξ", "u_η", "p", "T"],
               :Cond => [ "Dirichlet",
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
# second sub-domain

# θ = Cor[2]
# co = cos.(θ); sn = sin.(θ)
# left= Dict(:var => ["u_ξ", "u_η", "p", "T"],
#              :Cond => [ "Dirichlet",
#                         "Neumman",
#                         "Neumman",
#                         "Neumman"],
#              :value => [ 0,
#                          0,
#                          0,
#                          0]
#            )
# right = Dict(:var => ["u_ξ", "u_η", "p", "T"],
#                :Cond => [ "Dirichlet",
#                          "Dirichlet",
#                          "Dirichlet",
#                          "Dirichlet"],
#                :value => [ 692 .*co,
#                            -692 .*sn,
#                            101325,
#                            298]
#                )
#
# bot = Dict(:var => ["u_ξ", "u_η", "p", "T"],
#                :Cond => [ "Interfaz",
#                          "Interfaz",
#                          "Interfaz",
#                          "Interfaz"],
#                :value => [ 0,
#                            0,
#                            0,
#                            0]
#                )
#
# top = Dict(:var => ["u_ξ", "u_η", "p", "T"],
#                :Cond => [ "Neumman",
#                          "Neumman",
#                          "Neumman",
#                          "Neumman"],
#                :value => [ 0,
#                            0,
#                            0,
#                            0]
#                )
#
# push!(df, [left, right, bot, top])



return df
end
