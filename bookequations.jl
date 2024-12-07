function eq_3_15(M)
    base = 1 + ((gamma-1)/2)*M^2
    exponent = gamma/(1-gamma)
    return base^exponent #p/po
end

function eq_6_10(M,θ)#yeilds pressure ratio
    num = 2*gamma*(M^2)*sin(θ)^2
    den = gamma + 1
    term1 = num/den
    term2 = (gamma-1)/den
    return term1 - term2
end

function eq_6_17(M,θ)
    temp = ((gamma-1)/2)
    num1 = 1 + temp*(M^2)
    den1 = gamma*(M^2)*(sin(θ)^2) - temp
    num2 = (M^2)*(cos(θ)^2)
    den2 = 1 + temp*(M^2)*(sin(θ)^2)
    return sqrt((num1/den1) + (num2/den2))#M2
end

function eq_6_18(M,θ)
    num = (M^2)*(sin(θ)^2)-1
    den = ((gamma+1)/2)*M^2 - num
    tan_δ = cot(θ)*num/den
    return atan(tan_δ)#δ
end

function eq_6_24(M)
    term1 = 1/(gamma*M^2)
    term2 = ((gamma+1)/(4))*M^2
    term3 = (gamma+1)*(((gamma+1)/16)*M^4+((gamma-1)/2)*M^2+1)
    sinsquared = (term1)*(term2-1+sqrt(term3))#6.24
    return asin(sqrt(sinsquared)) #θmax
end

function eq_7_10(M)
    term1 = sqrt((gamma+1)/(gamma-1))
    term2 = sqrt(((gamma-1)/(gamma+1))*abs(M^2 - 1))#TODO: Check validity of abs()
    term3 = sqrt(abs(M^2 - 1))
    return term1*atan(term2) - atan(term3) #pm_angle
end