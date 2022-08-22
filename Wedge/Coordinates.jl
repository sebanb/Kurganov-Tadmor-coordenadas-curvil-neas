function Car2Cur(vecx, vecy, ξ, η, id, ndom)
    if id == "Polar" && ndom ==1
       θ = η+π
    elseif id == "Polar" && ndom == 3
      θ = (3/2) * π - deg2rad(15) + η
    elseif id == "Cartesian"
       θ = deg2rad(15)
    end
 vec_ξ = vecx*cos(θ) + vecy*sin(θ)
 vec_η = -vecx*sin(θ) + vecy*cos(θ)


 return vec_ξ, vec_η
end

function Cur2Car(vecξ, vecη, ξ, η, id)
   if id == "Polar"
      θ = π-η
   elseif id == "Cartesian"
       θ = deg2rad(15)
   end
   vec_x = vecξ*cos(θ) - vecη*sin(θ)
   vec_y = vecξ*sin(θ) + vecη*cos(θ)
   return vec_x, vec_y
end

function changcoor(ξ, η, id)
     if id == "Polar"
        x = ξ * cos(η)
        y = ξ * sin(η)

    elseif id == "Cartesian"
        θ = deg2rad(15)
        x = ξ * cos(θ) - η * sin(θ)
        y = ξ * sin(θ) + η * cos(θ)
    end
     return x, y
end
