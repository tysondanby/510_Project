function plotdetachvsM(Mlims,t,c,h)#default values
    Ms = collect(Mlims[1]:((Mlims[2]-Mlims[1])/(99)):Mlims[2])
    AoAmaxs = zeros(length(Ms))
    AoAmins = zeros(length(Ms))
    for i =1:1:length(Ms)
        AoAmaxs[i], AoAmins[i] = findsep(Ms[i],t,c,h)
    end
    detachplot = plot(Ms,[AoAmaxs AoAmins].*(180/pi), xlabel = "Mach Number",ylabel = "Angle of Attack (°)",  labels = ["Bottom Shock Detaches" "Top Shock Detaches"], legend = :outertopright)
    savefig(detachplot,"plots/detach/Detach vary M.png")
end

function plotdetachvst(M,tlims,c,h)#default values
    ts = collect(tlims[1]:((tlims[2]-tlims[1])/(99)):tlims[2])
    AoAmaxs = zeros(length(ts))
    AoAmins = zeros(length(ts))
    for i =1:1:length(ts)
        AoAmaxs[i], AoAmins[i] = findsep(M,ts[i],c,h)
    end
    detachplot = plot(ts,[AoAmaxs AoAmins].*(180/pi), xlabel = "Thickness to Chord Ratio (t/c)",ylabel = "Angle of Attack (°)",  labels = ["Bottom Shock Detaches" "Top Shock Detaches"], legend = :outertopright)
    savefig(detachplot,"plots/detach/Detach vary t.png")
end

function plotdetachvsc(M,t,clims,h)#default values
    cs = collect(clims[1]:((clims[2]-clims[1])/(99)):clims[2])
    AoAmaxs = zeros(length(cs))
    AoAmins = zeros(length(cs))
    for i =1:1:length(cs)
        AoAmaxs[i], AoAmins[i] = findsep(M,t,cs[i],h)
    end
    detachplot = plot(cs,[AoAmaxs AoAmins].*(180/pi), xlabel = "Maximum Thickness Location (d/c)",ylabel = "Angle of Attack (°)",  labels = ["Bottom Shock Detaches" "Top Shock Detaches"], legend = :outertopright)
    savefig(detachplot,"plots/detach/Detach vary d.png")
end

function plotdetachvsh(M,t,c,hlims)#default values
    hs = collect(hlims[1]:((hlims[2]-hlims[1])/(99)):hlims[2])
    AoAmaxs = zeros(length(hs))
    AoAmins = zeros(length(hs))
    for i =1:1:length(hs)
        AoAmaxs[i], AoAmins[i] = findsep(M,t,c,hs[i])
    end
    detachplot = plot(hs,[AoAmaxs AoAmins].*(180/pi), xlabel = "Camber Ratio (h/c)",ylabel = "Angle of Attack (°)",  labels = ["Bottom Shock Detaches" "Top Shock Detaches"], legend = :outertopright)
    savefig(detachplot,"plots/detach/Detach vary h.png")
end