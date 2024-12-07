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
    δmax = eq_6_18(M,θmax)
    δ1, δ2, δ3, δ4, A1,A2,A3,A4 = airfoilgeometry(t,c,h)
    AoAmin = δ1 - δmax
    AoAmax = δmax - δ2
    #println("AoA Bounds $AoAmin to $AoAmax")#TODO:debug
    return AoAmax, AoAmin
end

function weaksolution_6_18(δ,M) #TODO: use bracketing instead between 
    solvefunction(θ) = eq_6_18(M,θ) - δ
    return zeronewton(solvefunction,0.1) #TODO: 60 degrees shoud be a good guess
end

function pressures(Minf,P_inf,AoA,δ1, δ2, δ3, δ4)
    topturnangle = δ1 - AoA #TODO: Handle what happens when one of these is negative
    bottomturnangle = δ2 + AoA
    #=
    println("Top Turn Angle:")
    println(180*topturnangle/pi)
    println("Bottom Turn Angle:")
    println(180*bottomturnangle/pi)
    # =#
    P1 = 0.0
    P2 = 0.0
    M1 = 0.0
    M2 = 0.0
    topexpandangle = δ1 + δ3
    bottomexpandangle = δ2 + δ4
    #=
    println("Expansion Angles (top, bottom):")
    println(topexpandangle*180/pi)
    println(bottomexpandangle*180/pi)
    # =#
    if (topturnangle >= 0.0) && (bottomturnangle >= 0.0)
        topshockangle = weaksolution_6_18(topturnangle,Minf)
        bottomshockangle = weaksolution_6_18(bottomturnangle,Minf)
        #println("Case1")
        #=
        println("Top Shock Angle:")
        println(180*topshockangle/pi)
        println("Bottom Shock Angle:")
        println(180*bottomshockangle/pi)
        # =#
        P1 = eq_6_10(Minf,topshockangle)*P_inf#oblique pressure ratio
        P2 = eq_6_10(Minf,bottomshockangle)*P_inf
        M1 = eq_6_17(Minf,topshockangle) #obliquedownstreammach
        M2 = eq_6_17(Minf,bottomshockangle)
    elseif (topturnangle <= 0.0)
        bottomshockangle = weaksolution_6_18(bottomturnangle,Minf)
        #println("Case2")
        #=
        println("Bottom Shock Angle:")
        println(180*bottomshockangle/pi)
        # =#
        P2 = eq_6_10(Minf,bottomshockangle)*P_inf
        M2 = eq_6_17(Minf,bottomshockangle)
        PMinf = eq_7_10(Minf)#pm_angle
        PM1 = PMinf - topturnangle
        eq_7_10_0(Mach) = eq_7_10(Mach) - PM1
        M1 = zeronewton(eq_7_10_0,2)#machfrompm(PM1)
        Pinf_div_Po = eq_3_15(Minf)
        P1_div_Po = eq_3_15(M1)
        Po = 1/(Pinf_div_Po / P_inf)
        P1 = P1_div_Po * Po
    else
        topshockangle = weaksolution_6_18(topturnangle,Minf)
        #println("Case3")
        #=
        println("Top Shock Angle:")
        println(180*topshockangle/pi)
        # =#
        P1 = eq_6_10(Minf,topshockangle)*P_inf
        M1 = eq_6_17(Minf,topshockangle)
        PMinf = eq_7_10(Minf)#pm_angle
        PM2 = PMinf - bottomturnangle
        eq_7_10_00(Mach) = eq_7_10(Mach) - PM2
        M2 = zeronewton(eq_7_10_00,2)#machfrompm(PM1)
        Pinf_div_Po = eq_3_15(Minf)
        P2_div_Po = eq_3_15(M2)
        Po = 1/(Pinf_div_Po / P_inf)
        P2 = P2_div_Po * Po
    end

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
    #=
    println("Pressures:")
    println(P1)
    println(P2)
    println(P3)
    println(P4)
    # =#
    pressurevec = [P1,P2,P3,P4]
    Avec1 = A1*[sin(δ1-AoA),-cos(δ1-AoA)]
    Avec2 = A2*[sin(δ2+AoA),cos(δ2+AoA)]
    Avec3 = A3*[-sin(δ3+AoA),-cos(δ3+AoA)]
    Avec4 = A4*[-sin(δ4-AoA),cos(δ4-AoA)]
    Areanormals = [Avec1,Avec2,Avec3,Avec4]
    L = 0.0
    D = 0.0
    for i = 1:1:4
        D,L = [D,L] + Areanormals[i]*pressurevec[i]
    end
    return L,D
end

function analyzeairfoil(M,P,t,c,h)
    tolerance = 0.0
    range = 1.0-tolerance
    AoAmax, AoAmin = findsep(M,t,c,h)
    step = (AoAmax-AoAmin)*range/100
    AoAs = collect((range*AoAmin):step:(range*AoAmax))
    Ls = zeros(length(AoAs))
    Ds = zeros(length(AoAs))
    for i = 2:1:(length(AoAs)-1) #First and last points dont compute due to separation.
        L, D = airfoilliftdrag(M,P,AoAs[i],t,c,h)
        Ls[i] = L
        Ds[i] = D
    end
    return AoAs, Ls, Ds
end