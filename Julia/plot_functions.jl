
## Colors from https://commons.wikimedia.org/wiki/File:Cell_Cycle_2-2.svg
colI = "#f8951d" #orange
colG1 = "#94b6d2" #blue
colS = "#dc4040" #red
colG2 = "#698e6d" #green
colM = "#a67ba7" #purple
## Colort pallete https://personal.sron.nl/~pault/
blue1 = "#332288"
blue2 = "#88ccee"
blue3 = "#63bfeb"
green1 = "#44aa99"
green2 = "#117733"
green3 = "#82c3af"
yellow1 = "#999933"
yellow2 = "#ddcc77"
red1 = "#cc6677"
red2 = "#882255"
purple1 = "#aa4499"
gray1 = "#888888"

@pyimport mpl_toolkits.axes_grid1 as axgrid

function PlotSolution(sol, rad=false)
    clf()
    fig, ax1 = subplots()
    mylw = 2
    ax2 = ax1[:twinx]()
    ax1[:plot](sol.t, sol[3, :], color=red1, lw=mylw, label="Wee1") # Wee1
    # ax1[:plot](sol.t, sol[4, :], color=purple1, lw=mylw, label="Cdc25\$_P\$") # Cdc25_P
    ax1[:plot](sol.t, sol[6, :], color=gray1, lw=mylw, label="APC") # APC
    # ax1[:plot](sol.t, sol[2, :], color = "LimeGreen", lw=mylw, label="MPF\$_P\$") # MPF_P
    ax1[:plot](sol.t, sol[1, :], color=green2, lw=mylw, label="MPF") # MPF
    ax2[:plot](sol.t, sol[7, :], color="black", lw=mylw, label="Mass") #Mass 
    ax1[:plot]([sol.t[1], sol.t[end]], [-1, -1], color="black", lw=mylw, label="Mass") #Phantom Mass for lagend
    if rad
        ax1[:plot](sol.t, sol[8, :], "--", color=blue3, lw=mylw, label="Chk2", alpha=0.9) # Chk2
        ax1[:plot](sol.t, sol[9, :]/35., "--", color=green3, lw=mylw, label="ATM", alpha =0.6) # ATM
        ax1[:plot](sol.t, sol[10, :]/35., "--", color=blue2, lw=mylw, label="DSB", alpha =0.6) # DSB
        #=
        rad_start = sol.prob.f.params[41]
        rad_dur = sol.prob.f.params[42]
        rad_rate = sol.prob.f.params[43]
        rad_peak = [rad_start, sol(rad_start+rad_dur)[9]/35.0] + [0.0,0.02]
        str_rad = @sprintf("Radiation pulse\n\$D = %4.2f\$ Gy", rad_dur*rad_rate)
        arrowprops = Dict("arrowstyle"=>"fancy", "fc"=>"gray", "ec"=>"none")
        annotate(str_rad, xy=rad_peak, xytext=rad_peak+[0.0,0.2], ha="center", arrowprops=arrowprops)
        =#
    end
    ax1[:fill_between](sol.t, -1, -1, facecolor=colI, edgecolor="gray", alpha = 0.8, label="Interphase")
    ax1[:fill_between](sol.t, -1, -1, facecolor=colM, edgecolor="gray", alpha = 0.8, label="M-phase")
    leg = ax1[:legend](loc="upper left", bbox_to_anchor=(1.12, 1.0), fontsize=12, labelspacing=0.3, borderpad=0.4)
    leg[:set_zorder](6)
    ax1[:set_xlabel]("time (hr)", fontsize=14)
    ax1[:set_ylabel]("Concentration", fontsize=14)
    ax2[:set_ylabel]("Cell mass", fontsize=14)
    # str_Dose = @sprintf("\$D = %4.2f\$", TrapRule(map(DoseRate, t), t))
    # bbox_prop0 = Dict("boxstyle"=>"round", "fc"=>"white")
    # text(180.0, 1.25, str_Dose, size=12, color="DarkSlateGray", bbox=bbox_prop0)
    min = round(minimum(sol[1,:]), 1)
    max = round(maximum(sol[1,:]), 1)
    factor = (max - min)/100*4
    xlim(0.0, sol.t[end])
    ax1[:set_ylim](min - factor, max + factor)
    ax2[:set_ylim](0.0, sol.prob.p.K_Mass)
    
    # Cell cycle stages
    divider = axgrid.make_axes_locatable(ax1)
    ax0 = divider[:new_vertical](size="4%")
    fig[:add_axes](ax0)
    xx = sol.t
    yy = sol[1, :]
    yy2 = sol[6, :]
    # M-phase: MPF above threshold theta_M
    ixM = yy .> sol.prob.p.θ_M
    # Interphase MPF bellow threshold theta_M
    ixI = yy .< sol.prob.p.θ_M
    ax0[:fill_between](xx, 0, 1, where=ixI, facecolor=colI, edgecolor="gray", alpha = 0.8)
    ax0[:fill_between](xx, 0, 1, where=ixM, facecolor=colM, edgecolor="gray", alpha = 0.8)
    setp(ax0[:get_xticklabels](), visible=false)
    setp(ax0[:get_yticklabels](), visible=false)
    ax0[:xaxis][:set_ticks_position]("none")
    ax0[:yaxis][:set_ticks_position]("none")
    axis([0., sol.t[end], 0., 1.])
    sca(ax1)

    gcf()
end # PlotSol

function PlotBifurcation(curve::BifurcationCurve, vars)
    c = Dict(["S"=>"black", "N"=>"red", "U"=>"red", "LC"=>"blue"]) 
    ls = Dict(["S"=>"-", "N"=>"-.", "U"=>"--", "LC"=>"-"])
    lw = Dict(["S"=>2, "N"=>2, "U"=>2, "LC"=>3])
    len = length(curve.stab)
    x0 = 1
    x = curve.d[vars[1]]
    y = curve.d[vars[2]]
    for i in 1:length(curve.changes)
        xf = curve.changes[i]
        idx = x0:xf
        plot(x[idx], y[idx], c=c[curve.stab[x0+1]], ls=ls[curve.stab[x0+1]], lw=lw[curve.stab[x0+1]])
        x0 = xf
    end
    xf = len
    idx = x0:xf
    if typeof(y)<:Dict
        plot(x[idx], y[:min][idx], c=c[curve.stab[x0+1]], ls=ls[curve.stab[x0+1]], lw=lw[curve.stab[x0+1]])
        plot(x[idx], y[:max][idx], c=c[curve.stab[x0+1]], ls=ls[curve.stab[x0+1]], lw=lw[curve.stab[x0+1]])
    else
        plot(x[idx], y[idx], c=c[curve.stab[x0+1]], ls=ls[curve.stab[x0+1]], lw=lw[curve.stab[x0+1]])
    end
    points = curve.special_points
    for i in keys(points)
        for j in 1:length(points[i]) 
            if i in ["LP", "H", "BP"]
                plot(points[i][j][vars[1]], points[i][j][vars[2]], c="black", "o", ms=6)
                # text(points[i][j][vars[1]], points[i][j][vars[2]], i)
            end
        end
    end
end

