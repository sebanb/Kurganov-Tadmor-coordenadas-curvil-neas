function interfaz(ξ, η)
I = length(ξ)-1
J = length(η)-1
ξ_int = zeros(I)
η_int = zeros(J)
  for i in 1:I
      ξ_int[i] = ξ[i] + (ξ[i+1]-ξ[i])/2
  end
  for j=1:J
      η_int[j] = η[j] + (η[j+1]-η[j])/2
   end
return ξ_int, η_int
end

function metrics(ξ, η, id)
n = length(ξ)
m = length(η)
h_η = zeros(n, m); dhξ_dη = zeros(n, m);
h_ξ = zeros(n, m); dhη_dξ = zeros(n, m);
for i=1:n
   for j=1:m
    if id == "Polar"
        h_ξ[i, j] = 1
        h_η[i, j] = ξ[i]
        dhξ_dη[i, j] = 0
        dhη_dξ[i, j] = 1
    elseif id == "Parabolic"
       h_ξ[i, j] = sqrt(ξ[i]^2+η[j]^2)
       h_η[i, j] = sqrt(ξ[i]^2+η[j]^2)
       dhξ_dη[i, j] = η[j]/(sqrt(ξ[i]^2+η[j]^2))
       dhη_dξ[i, j] = ξ[i]/(sqrt(ξ[i]^2+η[j]^2))

    else
        h_ξ[i, j] = 1
        h_η[i, j] = 1
        dhξ_dη[i, j] = 0
        dhη_dξ[i, j] = 0
    end
   end
end
return h_ξ, h_η, dhξ_dη, dhη_dξ
end

function minmood(u, Δξ, Δη)
   I = size(u)[1]
   J = size(u)[2]
   u_x = zeros(I, J)
   u_y = zeros(I, J)
   for i=1:I
       for j=1:J
         if i == 1 # B.C in horizonal direction
           a1 = 0; b1 = 0; c1=0
        elseif i == I
           a1 = 0; b1 = 0; c1=0
         else
           a1 = (u[i,j]- u[i-1,j])/(Δξ)
           b1 = (u[i+1,j]- u[i-1,j])/(2*Δξ)
           c1 = (u[i+1,j]- u[i,j])/(Δξ)
         end
        if a1>0 && b1>0 && c1>0
           u_x[i, j] = min(a1, b1, c1)
        elseif a1<0 && b1<0 && c1<0
           u_x[i, j] = max(a1, b1, c1)
        else
           u_x[i,j] = 0
        end
        if j == 1 # B.C in vertical direction
           a2 = 0; b2 = 0; c2 = 0
        elseif j == J
           a2 = 0; b2 = 0; c2 = 0
         else
           a2 = (u[i,j]-u[i,j-1])/(Δη)
           b2 = (u[i,j+1]-u[i,j-1])/(2*Δη)
           c2 = (u[i,j+1]-u[i,j])/(Δη)
          end
          if a2>0 && b2>0 && c2>0
             u_y[i, j] = min(a2, b2, c2)
          elseif a2<0 && b2<0 && c2<0
             u_y[i, j] = max(a2, b2, c2)
          else
             u_y[i,j] = 0
          end
        end
   end
# Now we have the slopes in each cells. We'll find the uplus and uminus in
# interfaces
u_px = zeros(I-1, J)
u_mx = zeros(I-1, J)
u_py = zeros(I, J-1)
u_my = zeros(I, J-1)
 for i=1:I-1
   for j=1:J
      u_mx[i,j] = u[i,j] + u_x[i,j] * Δξ/2
      u_px[i,j] = u[i+1,j] - u_x[i+1,j] * Δξ/2
    end
  end
 for i in 1:I
   for j in 1:J-1
       u_my[i,j] = u[i,j] + u_y[i,j]*Δη/2
       u_py[i,j] = u[i,j+1] - u_y[i,j+1]*Δη/2
   end
 end
return u_px, u_mx, u_py, u_my
end

function soundspeed(γ, p, ρ)
   s = sqrt.(γ.*p./ρ)
   return s
end

function maxlocspeeds(u_p, u_m, s_p, s_m)
λ_1p = abs.(u_p .- s_p)
λ_2p = abs.(u_p)
λ_3p = abs.(u_p .+ s_p)
λ_1m = abs.(u_m .- s_m)
λ_2m = abs.(u_m)
λ_3m = abs.(u_m .+ s_m)
a_p = max.(max.(λ_1p, λ_2p, λ_3p), max.(λ_1m, λ_2m, λ_3m))
a_m = min.(min.(λ_1p, λ_2p, λ_3p), min.(λ_1m, λ_2m, λ_3m))
mx_fluxu = maximum([abs.(λ_3p) abs.(λ_3m) abs.(λ_1m) abs.(λ_1p)])
return a_p, a_m, mx_fluxu
end

function FluxF(u, v, ρ, p, E)
   F1_u = ρ.*u
   F2_u = ρ.*u.^2 .+ p
   F3_u = ρ.*u.*v
   F4_u = u.*(E .+ p)
   return F1_u, F2_u, F3_u, F4_u
end

function FluxG(u, v, ρ, p, E)
   G1_u = ρ.*v
   G2_u = ρ.*u.*v
   G3_u = ρ.*v.^2 .+ p
   G4_u = v.*(E .+ p)
   return  G1_u, G2_u, G3_u, G4_u
end

function source(u, v, ρ, p, h_ξ, h_η, dhξ_dη, dhη_dξ)

  met1 = dhξ_dη./(h_ξ .*h_η)
  met2 = dhη_dξ./(h_ξ .*h_η)

  S1 = met2.*(ρ.*v.^2 .+ p) .- ρ.*met1.*u.*v
  S2 = met1.*(ρ.*u.^2 .+ p) .- ρ.*met2.*u.*v

  return S1, S2
end




function Interfaz(Dom, Domn, Boundary, Boundary_n, var, Δ, Cor, id, var_d)
  n = size(var[Dom],1); m = size(var[Dom],2);
  if Boundary ==1
     row = 1; col= :; rowp = 2; colp = :; row_int1 = 1; col_int1 = :;
  elseif Boundary == 2
     row = n; col= :; rowp = n-1; colp = :; row_int1 = n-1; col_int1 = :;
  elseif Boundary == 3
     row = :; col= 1; rowp = :; colp = 2; row_int1 = :; col_int1 = 1;
  elseif Boundary == 4
     row = :; col= m; rowp = :; colp = m-1; row_int1 = :; col_int1 = m-1;
  end
  n_n = size(var[Domn],1); m_n = size(var[Domn],2)
  if Boundary_n == 1
            rown = 1; coln = :; row_intn = 1; col_intn = :;
  elseif Boundary_n == 2
            rown = n_n; coln = :; row_intn = n_n -1; col_intn = :;
  elseif Boundary_n == 3
            rown = :; coln = 1; row_intn = :; col_intn = 1;
  elseif Boundary_n == 4
            rown = :; coln = m_n; row_intn = :; col_intn = m_n-1;
  end
  Δξ_1  = Δ[2*Dom-1]; Δη_1 = Δ[2*Dom]
  Δξ_n = Δ[2*Domn-1]; Δη_n = Δ[2*Domn]
     h_ξ1, h_η1 = metrics(Cor[2*Dom-1], Cor[2*Dom], id[Dom,2])
     h_ξn, h_ηn = metrics(Cor[2*Domn-1], Cor[2*Domn], id[Domn,2])
     ξ_int1, η_int1 = interfaz(Cor[2*Dom-1], Cor[2*Dom])
     ξ_intn, η_intn = interfaz(Cor[2*Domn-1], Cor[2*Domn])
     if Boundary == 1 || Boundary == 2
        h_ξint1, h_ηint1 = metrics(ξ_int1, Cor[2*Dom] , id[Dom, 2])
     else
        h_ξint1, h_ηint1 = metrics(Cor[2*Dom-1], η_int1 , id[Dom, 2])
     end

     if Boundary_n == 1 || Boundary_n == 2
        h_ξintn, h_ηintn = metrics(ξ_intn, Cor[2*Domn], id[Domn, 2])
     else
         h_ξintn, h_ηintn = metrics(Cor[2*Domn-1], η_intn, id[Domn, 2])
     end
     V_1 = 0.5*Δξ_1*Δη_1 .*(
           h_ξint1[row_int1, col_int1] .*h_ηint1[row_int1, col_int1] .+
           h_ξ1[row, col] .*h_η1[row, col])
     V_n = 0.5*Δξ_n*Δη_n .*(
           h_ξintn[row_intn, col_intn] .*h_ηintn[row_intn, col_intn] .+
           h_ξn[rown, coln] .*h_ηn[rown, coln])
     if (Boundary ==1 && Boundary_n == 2) || (Boundary ==2 && Boundary_n ==1)
        var_n = var
     elseif  (Boundary ==3 && Boundary_n == 4) || (Boundary ==4 && Boundary_n == 3)
        var_n = var
     else
        var_n = var_d
     end
     var_2 = (V_1 .*var[Dom][row, col] +
                           V_n .*var_n[Domn][rown, coln])./(V_1 .+ V_n)

return var_2, row, col
end


function BC(Dom, Boundary, var, Condition, value)
   n = size(var[Dom],1); m = size(var[Dom],2);
   if Boundary ==1
      row = 1; col= :; rowp = 2; colp = :;
   elseif Boundary == 2
      row = n; col= :; rowp = n-1; colp = :;
   elseif Boundary == 3
      row = :; col= 1; rowp = :; colp = 2;
   elseif Boundary == 4
      row = :; col= m; rowp = :; colp = m-1;
   end
   if Condition == "Neumman"
      var[Dom][row, col] .= value .+ var[Dom][rowp, colp]
   elseif Condition == "Dirichlet"
      var[Dom][row, col] .= value
   end
   return var[Dom]
end

function quarterint(Cord, sign, dir, Δ)
   if dir == "x"
      Cordquart = Cord .+ sign .*Δ/4
   else
      Cordquart = Cord .+ sign .*Δ/4
   end
   return  Cordquart
end
