#points 1 and 4 are the front and back.
#t = thickness to chord ratio
#c = chord location of max thickness
#h = avg height of points 2 and 3
using Plots

global gamma = 1.4
global Tinf = 223.150
global Pinf = 26436.3
global rhoinf = 0.412707
global ainf = 299.463
global color1 = RGB(.9,0,0)#Red
global color2 = RGB(0,0,.7)#Blue
global color3 = RGB(.9,0,0)#RedRGB(.3,0,.2)#Dark port
global color4 = RGB(0,0,.7)#BlueRGB(0,0.4,0)#Deep green
global color5 = RGB(.9,0,0)#RedRGB(0.5,0.2,0) #burnt orange
global color6 = RGB(0,0,.7)#BlueRGB(0.5,0,0.3)#Maroon
global color7 = RGB(.9,0,0)#RedRGB(0,0.25,0.0)#Forest
global color8 = RGB(0,0,.7)#BlueRGB(0,0,0.7)#Blue

include("tools.jl")
include("bookequations.jl")
include("analysis.jl")
include("polarplots.jl")
include("detachplots.jl")

function plotscenario(M,Mlims,t,tlims,c,clims,h,hlims)#default conditions
    numberofseries = 7
    plotsvaryM(Mlims,t,c,h; numseries = numberofseries)
    plotsvaryt(M,tlims,c,h; numseries = numberofseries)
    plotsvaryc(M,t,clims,h; numseries = numberofseries)
    plotsvaryh(M,t,c,hlims; numseries = numberofseries)
    plotdetachvsM(Mlims,t,c,h)
    plotdetachvst(M,tlims,c,h)
    plotdetachvsc(M,t,clims,h)
    plotdetachvsh(M,t,c,hlims)
end

plotscenario(3,(1.5,4.5),0.15,(0,0.3),0.5,(.35,.65),0,(0.0,.075))