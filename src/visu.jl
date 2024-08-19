module Visu

#using ImpedancePredictionVIE
using ..ImpedancePredictionVIE #damit alles exportierte hier funktioniert
#import CompScienceMeshes: Mesh
using CompScienceMeshes
using Plots
using Plotly
using BEAST
using LinearAlgebra
using StaticArrays
#using Distributed

#import ..pointlist2xyzlist
#import ..realnodes
#import ..realvertices
#import ImpedancePredictionVIE.realvertices





function iplot(; size=(900,700))
    #plotlyjs()
    plotly()
    plt = Plots.plot(size = size, xlabel = "x", ylabel = "y", zlabel = "z")
    return plt
end

function iplot2(; kwargs...)
    plotly()
    plt = Plots.plot(; kwargs...)
    return plt
end


function points(pointlist, plt = iplot())
    Allx, Ally, Allz = pointlist2xyzlist(pointlist)
    scatter!(plt, Allx, Ally, Allz, markersize = 3, color = "magenta", label="")

    return plt
end


function mesh(m::Mesh, plt = iplot(); nodesize = 0.5, nodecolor = "blue", linewidth = 1.0, linecolor = "green")

    points = realvertices(m)
    Allx, Ally, Allz = pointlist2xyzlist(points)

    scatter!(plt, Allx, Ally, Allz, markersize = nodesize, color = nodecolor, label="")

    all_edges = skeleton(m,1).faces

    for edge in all_edges
        a = m.vertices[edge[1]]
        b = m.vertices[edge[2]]
        plot!(plt, [a[1], b[1]], [a[2], b[2]], [a[3], b[3]], lcolor=linecolor, linewidth = linewidth, label="")
    end

    return plt
end

function meshp(m::Mesh, plt = iplot() )# = nothing)

    points = realvertices(m)
    Allx, Ally, Allz = pointlist2xyzlist(points)

    scatter!(plt, Allx, Ally, Allz, markersize = 1, color = "blue", label="")

    return plt
end


function fnspos(X::BEAST.Space, plt = iplot())

    active_positions = Vector{SVector{3, Float64}}()
    
    for (i,f) in enumerate(X.fns)
        !isempty(f) && push!(active_positions, X.pos[i]) 
    end

    Allx, Ally, Allz = pointlist2xyzlist(active_positions)
    scatter!(plt, Allx, Ally, Allz, markersize = 2.2, color ="red", label="")
    
    return plt
end


function add1(plt, simplex0, refsp, scale, locfnsnr, locsign) #adds one cell field...

    bary_list = [[u, v, w] for u in 0.1:0.3:1.0 for v in 0.1:0.3:1.0-u for w in 0.1:0.3:1.0-u-v]
    #bary_list = [[u, v, w] for u in 0.0:0.3:1.0 for v in 0.0:0.26:1.0-u for w in 0.0:0.26:1.0-u-v]

    # b1 = 0.25
    # b2 = 0.5
    # b3 = 0.75
    # bary_list = [
    # [b1, b1, b1],
    # [b1, b1, b2],
    # [b1, b1, b3],
    # [b1, b2, b1],
    # [b1, b2, b2],
    # [b1, b3, b1],
    # [b2, b1, b1],
    # [b2, b1, b2],
    # [b2, b2, b1],
    # [b3, b1, b1],
    # ]

    for bary in bary_list
        @assert norm(bary) <= 1.0
        n = neighborhood(simplex0, bary) #CSM-MP
        p = n.cart
        v = refsp(n)[locfnsnr].value*locsign  #<------- achtung tet hat 4 flächen muss man später spezifizieren.... erstmal die erste
        @assert length(v) == 3
        draw_arrow!(plt, p, v; scale = scale, arrcolor = "red", arrwidth = 2.5)
    end


    return plt
end


function simplex(plt, simplex0::CompScienceMeshes.Simplex; linewidth = 1.0, color = "red") #adds one cell field...

    if length(simplex0.vertices) == 4
        a = simplex0.vertices[1]
        b = simplex0.vertices[2]
        c = simplex0.vertices[3]
        d = simplex0.vertices[4]

        plt = plot!(plt,[a[1], b[1]], [a[2], b[2]], [a[3], b[3]], linewidth=linewidth, color=color, label="")
        plt = plot!(plt,[a[1], c[1]], [a[2], c[2]], [a[3], c[3]], linewidth=linewidth, color=color, label="")
        plt = plot!(plt,[a[1], d[1]], [a[2], d[2]], [a[3], d[3]], linewidth=linewidth, color=color, label="")


        plt = plot!(plt,[b[1], c[1]], [b[2], c[2]], [b[3], c[3]], linewidth=linewidth, color=color, label="")
        plt = plot!(plt,[b[1], d[1]], [b[2], d[2]], [b[3], d[3]], linewidth=linewidth, color=color, label="")

        plt = plot!(plt,[c[1], d[1]], [c[2], d[2]], [c[3], d[3]], linewidth=linewidth, color=color, label="")

    elseif length(simplex0.vertices) == 3
        a = simplex0.vertices[1]
        b = simplex0.vertices[2]
        c = simplex0.vertices[3]
    
        plt = plot!(plt,[a[1], b[1]], [a[2], b[2]], [a[3], b[3]], linewidth=linewidth, color=color, label="")
        plt = plot!(plt,[a[1], c[1]], [a[2], c[2]], [a[3], c[3]], linewidth=linewidth, color=color, label="")
    
        plt = plot!(plt,[b[1], c[1]], [b[2], c[2]], [b[3], c[3]], linewidth=linewidth, color=color, label="")
    else
        error()
    end

    return plt
end

function simplex(plt, fourpoints) #adds one cell field...
    @assert length(fourpoints) == 4
    a = fourpoints[1]
    b = fourpoints[2]
    c = fourpoints[3]
    d = fourpoints[4]

    plt = plot!(plt,[a[1], b[1]], [a[2], b[2]], [a[3], b[3]], color="red", label="")
    plt = plot!(plt,[a[1], c[1]], [a[2], c[2]], [a[3], c[3]], color="red", label="")
    plt = plot!(plt,[a[1], d[1]], [a[2], d[2]], [a[3], d[3]], color="red", label="")


    plt = plot!(plt,[b[1], c[1]], [b[2], c[2]], [b[3], c[3]], color="red", label="")
    plt = plot!(plt,[b[1], d[1]], [b[2], d[2]], [b[3], d[3]], color="red", label="")

    plt = plot!(plt,[c[1], d[1]], [c[2], d[2]], [c[3], d[3]], color="red", label="")

    return plt
end




function draw_arrow!(plt, pnt, v; scale = 1.0, arrcolor = "red", arrwidth = 1)
    # Funktion zum Zeichnen eines Pfeils von Punkt a nach Punkt b

    a = pnt
    b = pnt + v*scale

    plot!(plt, [a[1], b[1]], [a[2], b[2]], [a[3], b[3]], color=arrcolor, label="", linewidth=arrwidth)

    # Berechnen Sie die Punkte für die Pfeilspitze
    
    r = rand(3)
    k1 = (a - b)/norm(a - b) # k1 und k2 sind nicht orthogonale Einheitsvektoren
    k2 = (r - b)/norm(r - b)

    # orthogonalisiere k2 bzgl. k1
    k2 = k2 - dot(k1,k2)*k1
    #@assert dot(k1,k2) < 1e-8
    nrm = norm(a - b)
    d = b + k1*nrm/5 + k2*nrm/7    # Achtung k1 ist Rückwärts!
    e = b + k1*nrm/5 - k2*nrm/7

    # Zeichne Pfeilspitze
    plot!(plt, [b[1], d[1]], [b[2], d[2]], [b[3], d[3]], color=arrcolor, label="", linewidth=arrwidth)
    plot!(plt, [b[1], e[1]], [b[2], e[2]], [b[3], e[3]], color=arrcolor, label="", linewidth=arrwidth)
    return plt
end



function fieldplot(pointlist, vectorlist, scale, plt = iplot())
    #avgnorm = sum(norm.(vectorlist))/length(vectorlist)
    #scale = scale*1/avgnorm  ja....


    scale = scale*1/maximum(norm.(vectorlist))
    #@show scale
    @assert length(pointlist) == length(vectorlist)

    for i = 1:1:length(pointlist)
    draw_arrow!(plt, pointlist[i], vectorlist[i], scale = scale, arrcolor = "red", arrwidth = 2)
    end

    return plt
end


end #module