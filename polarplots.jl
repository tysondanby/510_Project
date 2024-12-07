function plotsvaryM(Mlims,t,c,h; numseries = 5, directory = "plots/polars/")
    Ms = collect(Mlims[1]:((Mlims[2]-Mlims[1])/(numseries - 1)):Mlims[2])
    AoAseries = []#TODO: declare variable type
    Lseries = []
    Dseries = []
    L_Dseries = []
    Clseries = []
    Cdseries = []
    seriescolors = Matrix{RGB}(zeros(1,numseries))
    serieslabels = Matrix{String}(undef,1,numseries)
    for i = 1:1:numseries
        M = Ms[i]
        x = (i-1)/(numseries-1)
        vinf = M*ainf
        AoAs, Ls, Ds = analyzeairfoil(M,Pinf,t,c,h)
        L_Ds = Ls ./ Ds
        Cls = Ls / (0.5*rhoinf*vinf^2) #A = 1
        Cds = Ds / (0.5*rhoinf*vinf^2)
        push!(AoAseries,AoAs[2:end-1]*180/pi)
        push!(Lseries,Ls[2:end-1])
        push!(Dseries,Ds[2:end-1])
        push!(L_Dseries,L_Ds[2:end-1])
        push!(Clseries,Cls[2:end-1])
        push!(Cdseries,Cds[2:end-1])
        seriescolors[i] = (1-x)*color1 + x*color2
        serieslabels[i] = "M = $(round(M,digits=2))"
    end
    Lplot = plot(AoAseries,Lseries, xlabel = "Angle of Attack (°)",ylabel = "Lift (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Dplot = plot(AoAseries,Dseries, xlabel = "Angle of Attack (°)",ylabel = "Drag (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    L_Dplot = plot(AoAseries,L_Dseries, xlabel = "Angle of Attack (°)",ylabel = "Lift / Drag", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Clplot = plot(AoAseries,Clseries, xlabel = "Angle of Attack (°)",ylabel = "Lift Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Cdplot = plot(AoAseries,Cdseries, xlabel = "Angle of Attack (°)",ylabel = "Drag Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    savefig(Lplot,directory*"L vary M.png")
    savefig(Dplot,directory*"D vary M.png")
    savefig(L_Dplot,directory*"L_D vary M.png")
    savefig(Clplot,directory*"Cl vary M.png")
    savefig(Cdplot,directory*"Cd vary M.png")
end

function plotsvaryt(M,tlims,c,h; numseries = 5, directory = "plots/polars/") #changes between plottingfunc
    ts = collect(tlims[1]:((tlims[2]-tlims[1])/(numseries - 1)):tlims[2])#changes between plottingfunc
    AoAseries = []#TODO: declare variable type
    Lseries = []
    Dseries = []
    L_Dseries = []
    Clseries = []
    Cdseries = []
    seriescolors = Matrix{RGB}(zeros(1,numseries))
    serieslabels = Matrix{String}(undef,1,numseries)
    vinf = M*ainf
    for i = 1:1:numseries
        t = ts[i] #changes between plottingfunc
        println("Plotting t = $t")
        x = (i-1)/(numseries-1)
        AoAs, Ls, Ds = analyzeairfoil(M,Pinf,t,c,h)
        L_Ds = Ls ./ Ds
        Cls = Ls / (0.5*rhoinf*vinf^2) #A = 1
        Cds = Ds / (0.5*rhoinf*vinf^2)
        push!(AoAseries,AoAs[2:end-1]*180/pi)
        push!(Lseries,Ls[2:end-1])
        push!(Dseries,Ds[2:end-1])
        push!(L_Dseries,L_Ds[2:end-1])
        push!(Clseries,Cls[2:end-1])
        push!(Cdseries,Cds[2:end-1])
        seriescolors[i] = (1-x)*color3 + x*color4 #changes between plottingfunc
        serieslabels[i] = "t/c = $(round(t,digits=2))" #changes between plottingfunc
    end
    Lplot = plot(AoAseries,Lseries, xlabel = "Angle of Attack (°)",ylabel = "Lift (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Dplot = plot(AoAseries,Dseries, xlabel = "Angle of Attack (°)",ylabel = "Drag (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    L_Dplot = plot(AoAseries,L_Dseries, xlabel = "Angle of Attack (°)",ylabel = "Lift / Drag", c = seriescolors, labels = serieslabels,ylims = (-25,25), legend = :outertopright)
    Clplot = plot(AoAseries,Clseries, xlabel = "Angle of Attack (°)",ylabel = "Lift Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Cdplot = plot(AoAseries,Cdseries, xlabel = "Angle of Attack (°)",ylabel = "Drag Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    endofstring = " vary t.png" #changes between plottingfunc
    savefig(Lplot,directory*"L"*endofstring)
    savefig(Dplot,directory*"D"*endofstring)
    savefig(L_Dplot,directory*"L_D"*endofstring)
    savefig(Clplot,directory*"Cl"*endofstring)
    savefig(Cdplot,directory*"Cd"*endofstring)
end

function plotsvaryc(M,t,clims,h; numseries = 5, directory = "plots/polars/") #changes between plottingfunc
    cs = collect(clims[1]:((clims[2]-clims[1])/(numseries - 1)):clims[2])#changes between plottingfunc
    AoAseries = []#TODO: declare variable type
    Lseries = []
    Dseries = []
    L_Dseries = []
    Clseries = []
    Cdseries = []
    seriescolors = Matrix{RGB}(zeros(1,numseries))
    serieslabels = Matrix{String}(undef,1,numseries)
    vinf = M*ainf
    for i = 1:1:numseries
        c = cs[i] #changes between plottingfunc
        println("Plotting c = $c")
        x = (i-1)/(numseries-1)
        AoAs, Ls, Ds = analyzeairfoil(M,Pinf,t,c,h)
        L_Ds = Ls ./ Ds
        Cls = Ls / (0.5*rhoinf*vinf^2) #A = 1
        Cds = Ds / (0.5*rhoinf*vinf^2)
        push!(AoAseries,AoAs[2:end-1]*180/pi)
        push!(Lseries,Ls[2:end-1])
        push!(Dseries,Ds[2:end-1])
        push!(L_Dseries,L_Ds[2:end-1])
        push!(Clseries,Cls[2:end-1])
        push!(Cdseries,Cds[2:end-1])
        seriescolors[i] = (1-x)*color5 + x*color6 #changes between plottingfunc
        serieslabels[i] = "d/c = $(round(c,digits=2))" #changes between plottingfunc
    end
    Lplot = plot(AoAseries,Lseries, xlabel = "Angle of Attack (°)",ylabel = "Lift (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Dplot = plot(AoAseries,Dseries, xlabel = "Angle of Attack (°)",ylabel = "Drag (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    L_Dplot = plot(AoAseries,L_Dseries, xlabel = "Angle of Attack (°)",ylabel = "Lift / Drag", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Clplot = plot(AoAseries,Clseries, xlabel = "Angle of Attack (°)",ylabel = "Lift Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Cdplot = plot(AoAseries,Cdseries, xlabel = "Angle of Attack (°)",ylabel = "Drag Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    endofstring = " vary d.png" #changes between plottingfunc
    savefig(Lplot,directory*"L"*endofstring)
    savefig(Dplot,directory*"D"*endofstring)
    savefig(L_Dplot,directory*"L_D"*endofstring)
    savefig(Clplot,directory*"Cl"*endofstring)
    savefig(Cdplot,directory*"Cd"*endofstring)
end

function plotsvaryh(M,t,c,hlims; numseries = 5, directory = "plots/polars/") #changes between plottingfunc
    hs = collect(hlims[1]:((hlims[2]-hlims[1])/(numseries - 1)):hlims[2])#changes between plottingfunc
    AoAseries = []#TODO: declare variable type
    Lseries = []
    Dseries = []
    L_Dseries = []
    Clseries = []
    Cdseries = []
    seriescolors = Matrix{RGB}(zeros(1,numseries))
    serieslabels = Matrix{String}(undef,1,numseries)
    vinf = M*ainf
    for i = 1:1:numseries
        h = hs[i] #changes between plottingfunc
        println("Plotting h = $h")
        x = (i-1)/(numseries-1)
        AoAs, Ls, Ds = analyzeairfoil(M,Pinf,t,c,h)
        L_Ds = Ls ./ Ds
        Cls = Ls / (0.5*rhoinf*vinf^2) #A = 1
        Cds = Ds / (0.5*rhoinf*vinf^2)
        push!(AoAseries,AoAs[2:end-1]*180/pi)
        push!(Lseries,Ls[2:end-1])
        push!(Dseries,Ds[2:end-1])
        push!(L_Dseries,L_Ds[2:end-1])
        push!(Clseries,Cls[2:end-1])
        push!(Cdseries,Cds[2:end-1])
        seriescolors[i] = (1-x)*color7 + x*color8 #changes between plottingfunc
        serieslabels[i] = "h/c = $(round(h,digits=3))" #changes between plottingfunc
    end
    Lplot = plot(AoAseries,Lseries, xlabel = "Angle of Attack (°)",ylabel = "Lift (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Dplot = plot(AoAseries,Dseries, xlabel = "Angle of Attack (°)",ylabel = "Drag (N)", c = seriescolors, labels = serieslabels, legend = :outertopright)
    L_Dplot = plot(AoAseries,L_Dseries, xlabel = "Angle of Attack (°)",ylabel = "Lift / Drag", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Clplot = plot(AoAseries,Clseries, xlabel = "Angle of Attack (°)",ylabel = "Lift Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    Cdplot = plot(AoAseries,Cdseries, xlabel = "Angle of Attack (°)",ylabel = "Drag Coefficient", c = seriescolors, labels = serieslabels, legend = :outertopright)
    endofstring = " vary h.png" #changes between plottingfunc
    savefig(Lplot,directory*"L"*endofstring)
    savefig(Dplot,directory*"D"*endofstring)
    savefig(L_Dplot,directory*"L_D"*endofstring)
    savefig(Clplot,directory*"Cl"*endofstring)
    savefig(Cdplot,directory*"Cd"*endofstring)
end