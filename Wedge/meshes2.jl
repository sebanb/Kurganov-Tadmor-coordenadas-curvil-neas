# This file is for build the mesh

function Mesh()

theta = 15 # Wedge angle
L = 1 # wedge long

## -----------------------
# 1st sub-domain(Polar)
## -----------------------
r_1 = 0.5 # Radius
theta_1 = (90-theta)*pi/180

nξ_1 = 120 # Number of elements which it's splited the radius
nη_1  = 80 # Number of elements which it's splited the angle

ξ_1 = LinRange(0, r_1, nξ_1)
η_1 = LinRange(0, theta_1, nη_1 )

Δξ_1 = r_1/(nξ_1 - 1)
Δη_1 = (theta_1)/(nη_1  - 1)
bc1 = [1 "Path" 0 0; 2 "Path" 0 0; 3 "Path" 0 0; 4 "Interfaz" 2 1]



## -----------------------
# 2nd sub-domain(Polar)
## -----------------------

diagwedge = L/(cos(deg2rad(theta)))
nξ_2 = 130
nη_2 = nξ_1

ξ_2 = LinRange(0, diagwedge, nξ_2)
η_2 = LinRange(0, r_1, nη_2)

Δξ_2 = diagwedge/(nξ_2  - 1)
Δη_2 = (r_1)/(nη_2   - 1)
bc2 = [1 "Interfaz" 1 4; 2 "Interfaz" 3 3 ; 3 "Path" 0 0; 4 "Path" 0 0]

## -----------------------
# 3rd sub-domain(Polar)
## -----------------------
nξ_3 = nη_2
nη_3 = 10
theta_3 = theta*π/180
ξ_3 = LinRange(0, r_1, nξ_3)
η_3 = LinRange(0, theta_3 , nη_3)

Δξ_3 = r_1/(nξ_3  - 1)
Δη_3 = (theta_3)/(nη_3   - 1)
bc3 = [1 "Path"  0 0; 2 "Path" 0 0; 3 "Interfaz" 2 2; 4 "Path" 0 0]
## ----------------------------------------------------------------------------

Dom = [1 "Polar"; 2 "Cartesian"; 3 "Polar"]

Cor = [ξ_1, η_1, ξ_2, η_2, ξ_3, η_3]

Δ = [Δξ_1, Δη_1, Δξ_2, Δη_2, Δξ_3, Δη_3]

BC = [bc1, bc2, bc3]

return Dom, Cor, Δ, BC
end
