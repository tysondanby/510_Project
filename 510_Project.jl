#points 1 and 4 are the front and back.
#t = thickness to chord ratio
#c = chord location of max thickness
#h = avg height of points 2 and 3
global gamma = 1.4

function newtonitter(x,f)
    dx = 1e-7
    dfdx = (f(x+dx)-f(x-dx))/(2*dx)
    return x - f(x)/dfdx
end

function zeronewton(f,xi)
    xnew = copy(xi)
    xold = 2*xnew + 1
    while abs((xnew-xold)/xold) > 1e-3
        #println(xold)
        #println(xnew)
        xold = copy(xnew)
        xnew = newtonitter(xnew,f)
    end
    return xnew
end

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
    term2 = sqrt(((gamma-1)/(gamma+1))*(M^2 - 1))
    term3 = sqrt(M^2 - 1)
    return term1*atan(term2) - atan(term3) #pm_angle
end

function airfoilgeometry(t,c,h)
    δ1 = atan((h+0.5*t)/c)
    δ2 = atan((0.5*t-h)/c)
    δ3 = atan((h+0.5*t)/(1-c))
    δ4 = atan((0.5*t-h)/(1-c))
    A1 = sqrt((0.5*t+h)^2+c^2)
    A2 = sqrt((0.5*t-h)^2+c^2)
    A3 = sqrt((0.5*t+h)^2+(1-c)^2)
    A4 = sqrt((0.5*t-h)^2+(1-c)^2)
    return δ1, δ2, δ3, δ4, A1,A2,A3,A4
end

function findsep(M,t,c,h)
    θmax = eq_6_24(M)
    δmax = eq_6_18(M,θ)
    δ1, δ2, δ3, δ4, A1,A2,A3,A4 = airfoilgeometry(t,c,h)
    AoAmin = δ1 - δmax
    AoAmax = δmax - δ2
    return AoAmax, AoAmin
end

function weaksolution_6_18(δ,M) #TODO: use bracketing instead between 
    solvefunction(θ) = eq_6_18(M,θ) - δ
    return zeronewton(solvefunction,0.1) #TODO: 60 degrees shoud be a good guess
end

function pressures(Minf,Pinf,AoA,δ1, δ2, δ3, δ4)
    topturnangle = δ1 - AoA #TODO: assert that these truning angles are positive.
    println("Top Turn Angle:")
    println(180*topturnangle/pi)
    bottomturnangle = δ1 + AoA
    println("Bottom Turn Angle:")
    println(180*bottomturnangle/pi)
    topexpandangle = δ1 + δ3
    bottomexpandangle = δ2 + δ4
    println("Expansion Angles (top, bottom):")
    println(topexpandangle*180/pi)
    println(bottomexpandangle*180/pi)
    topshockangle = weaksolution_6_18(topturnangle,Minf)
    bottomshockangle = weaksolution_6_18(bottomturnangle,Minf)
    println("Top Shock Angle:")
    println(180*topshockangle/pi)
    println("Bottom Shock Angle:")
    println(180*bottomshockangle/pi)
    P1 = eq_6_10(Minf,topshockangle)*Pinf#oblique pressure ratio
    P2 = eq_6_10(Minf,bottomshockangle)*Pinf
    M1 = eq_6_17(Minf,topshockangle) #obliquedownstreammach
    M2 = eq_6_17(Minf,bottomshockangle)
    PM1 = eq_7_10(M1)#pm_angle
    PM2 = eq_7_10(M2)
    PM3 = PM1 + topexpandangle
    PM4 = PM2 + bottomexpandangle
    eq_7_10_1(Mach) = eq_7_10(Mach) - PM3
    eq_7_10_2(Mach) = eq_7_10(Mach) - PM4
    M3 = zeronewton(eq_7_10_1,2)#machfrompm(PM3)
    M4 = zeronewton(eq_7_10_2,2)
    P1_div_Po = eq_3_15(M1)
    P2_div_Po = eq_3_15(M2)
    P3_div_Po = eq_3_15(M3)
    P4_div_Po = eq_3_15(M4)
    Potop = 1/(P1_div_Po / P1)
    Pobottom = 1/(P2_div_Po / P2)
    P3 = P3_div_Po * Potop
    P4 = P4_div_Po * Pobottom
    return P1,P2,P3,P4
end

function airfoilliftdrag(M,P,AoA,t,c,h)#TODO
    δ1, δ2, δ3, δ4, A1,A2,A3,A4 = airfoilgeometry(t,c,h)
    P1,P2,P3,P4 = pressures(M,P,AoA,δ1, δ2, δ3, δ4)
    println("Pressures:")
    println(P1)
    println(P2)
    println(P3)
    println(P4)
end

function analyzeairfoil(M,t,c,h)
    AoAmax, AoAmin = findsep(t,c,h)
    step = (AoAmax-AoAmin)/100
    AoAs = collect(AoAmin::AoAmax)
    Ls = zeros(length(AoAs))
    Ds = zeros(length(AoAs))
    @. Ls, Ds = airfoilliftdrag(AoAs,t,c,h)
    return AoAs, Ls, Ds
end