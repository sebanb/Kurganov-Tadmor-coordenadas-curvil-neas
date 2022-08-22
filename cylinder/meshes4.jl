
function Mesh()
#-------------------------------------
r = 1
rmesh = 8*r
n1r = 75
n1a = 80

ξ_1 = LinRange(r, rmesh, n1r)
η_1 = LinRange(0, π/2, n1a)

Δξ_1 = (rmesh - r)/(n1r-1)
Δη_1 = (π/2)/(n1a-1)
#-------------------------------------
# ξ_2 = LinRange(0, rmesh, n1r)
# η_2 = LinRange(π/2, π, n1a)
#
# Δξ_2 = rmesh/(n1r-1)
# Δη_2 = (π/2)/(n1a-1)



Dom = [1 "Polar"]#; 2 "Polar"]
Cor = [ξ_1,  η_1]#, ξ_2, η_2]
Δ = [Δξ_1,  Δη_1]#, Δξ_2, Δη_2]
bc = [1 "Path" 0 0; 2 "Path" 0 0; 3 "Path" 0 0; 4 "Path" 0 0]
BC = [bc]
return Dom, Cor, Δ, BC
end
