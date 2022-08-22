
function Mesh()
#-------------------------------------
vcte_wall = 0.77
vcte_sup = sqrt(5)
ucte_ini = 0
ucte_sup = sqrt(6.7449)
nu = 90
nv = 90

ξ_1 = LinRange(ucte_ini, ucte_sup, nu)
η_1 = LinRange(vcte_wall, vcte_sup, nv)

Δξ_1 = (ucte_sup - ucte_ini)/(nu-1)
Δη_1 = (vcte_sup - vcte_wall)/(nv-1)
#-------------------------------------
# ξ_2 = LinRange(0, rmesh, n1r)
# η_2 = LinRange(π/2, π, n1a)
#
# Δξ_2 = rmesh/(n1r-1)
# Δη_2 = (π/2)/(n1a-1)



Dom = [1 "Parabolic"]#; 2 "Polar"]
Cor = [ξ_1,  η_1]#, ξ_2, η_2]
Δ = [Δξ_1,  Δη_1]#, Δξ_2, Δη_2]
bc = [1 "Path" 0 0; 2 "Path" 0 0; 3 "Path" 0 0; 4 "Path" 0 0]
BC = [bc]
return Dom, Cor, Δ, BC
end
