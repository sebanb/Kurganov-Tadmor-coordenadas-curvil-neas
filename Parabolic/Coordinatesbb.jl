function Car2Cur(vecx, vecy, ξ, η, id)
    if id == "Polar"
       θ = η+π
       vec_ξ = vecx*cos(θ) + vecy*sin(θ)
       vec_η = -vecx*sin(θ) + vecy*cos(θ)
    elseif id == "Cartesian"
       θ = deg2rad(15)
       vec_ξ = vecx*cos(θ) + vecy*sin(θ)
       vec_η = -vecx*sin(θ) + vecy*cos(θ)
    elseif id == "Parabolic"
       θ = deg2rad(15)
       vec_ξ = ξ*vecx/(sqrt(ξ^2+η^2)) +  η*vecy/(sqrt(ξ^2+η^2))
       vec_η = -η*vecx/(sqrt(ξ^2+η^2)) +  ξ*vecy/(sqrt(ξ^2+η^2))
    end
 return vec_ξ, vec_η
end

function Cur2Car(vecξ, vecη, ξ, η, id)
   if id == "Polar"
      θ = π-η
      vec_x = vecξ*cos(θ) - vecη*sin(θ)
      vec_y = vecξ*sin(θ) + vecη*cos(θ)
   elseif id == "Cartesian"
       θ = deg2rad(15)
       vec_x = vecξ*cos(θ) - vecη*sin(θ)
       vec_y = vecξ*sin(θ) + vecη*cos(θ)
   elseif id == "Parabolic"
        θ = deg2rad(15)
        vec_x = vecξ*cos(θ) - vecη*sin(θ)
        vec_y = vecξ*sin(θ) + vecη*cos(θ)
   end

   return vec_x, vec_y
end

function changcoor(ξ, η, id)
     if id == "Polar"
        x = ξ * cos(η)
        y = ξ * sin(η)

    elseif id == "Parabolic"
        x = 0.5*(ξ^2 - η^2) 
        y = ξ*η
    end
     return x, y
end
