function visu(m::Mesh)

    # some meshes have "inactive" vertices, that do not contribute to the faces 
    # remove them!

    # m.faces
    # active_nodes=[] #index list

    # for (i, vi) in enumerate(m.vertices)
    #     breakfaceloop = false
    #     for (j, fj) in enumerate(m.faces)
    #         for n in fj
    #             if n==i
    #                 push!(active_nodes, i)
    #                 breakfaceloop = true
    #                 break
    #             end
    #         end
    #         breakfaceloop && break
    #     end
    # end

    plotly()
    AllX=[]
    AllY=[]
    AllZ=[]
    #active_vertices = m.vertices[active_nodes] 
    active_nodes = []
    for n in skeleton(m,0).faces
        push!(active_nodes,n[1])
    end 
    active_vertices = m.vertices[active_nodes]
    for p in active_vertices
        push!(AllX, p[1])
        push!(AllY, p[2])
        push!(AllZ, p[3])
    end
    plt=scatter(AllX, AllY, AllZ, markersize = 1, size=(850,850))
    plt=plot!(plt, xlabel="x", ylabel="y", zlabel="z")

    # Plot edges, igore double plot! NE IGNORIEREN KÖNNTE MASSIVST GRAFIKLEITUNG FORDERN
    # Folglich darf eine neue kante nur ergänzt werden wenn sie nicht schon ergänzt wurde


    all_edges = skeleton(m,1).faces
    for edge in all_edges
        a = m.vertices[edge[1]]
        b = m.vertices[edge[2]]
        plt = plot!(plt,[a[1], b[1]], [a[2], b[2]], [a[3], b[3]], color="green", label="")
    end

    # if length(m.faces[1]) == 4
    #     for face in m.faces
    #         a = m.vertices[face[1]]
    #         b = m.vertices[face[2]]
    #         c = m.vertices[face[3]]
    #         d = m.vertices[face[4]]

    #         plt = plot!(plt,[a[1], b[1]], [a[2], b[2]], [a[3], b[3]], color="green", label="")
    #         plt = plot!(plt,[a[1], c[1]], [a[2], c[2]], [a[3], c[3]], color="green", label="")
    #         plt = plot!(plt,[a[1], d[1]], [a[2], d[2]], [a[3], d[3]], color="green", label="")


    #         plt = plot!(plt,[b[1], c[1]], [b[2], c[2]], [b[3], c[3]], color="green", label="")
    #         plt = plot!(plt,[b[1], d[1]], [b[2], d[2]], [b[3], d[3]], color="green", label="")

    #         plt = plot!(plt,[c[1], d[1]], [c[2], d[2]], [c[3], d[3]], color="green", label="")
    #     end
    # elseif length(m.faces[1]) == 3
    #     for face in m.faces
    #         a = m.vertices[face[1]]
    #         b = m.vertices[face[2]]
    #         c = m.vertices[face[3]]

    #         plt = plot!(plt,[a[1], b[1]], [a[2], b[2]], [a[3], b[3]], color="green", label="")
    #         plt = plot!(plt,[a[1], c[1]], [a[2], c[2]], [a[3], c[3]], color="green", label="")

    #         plt = plot!(plt,[b[1], c[1]], [b[2], c[2]], [b[3], c[3]], color="green", label="")

    #     end
    # else 
    #     error()
    # end


    return plt
end



# function visu(ANDERER TYP)

#end