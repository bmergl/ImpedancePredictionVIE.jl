


a = [i^2 for i in 1:100000000]

b = [(i+10)^2 for i in 1:10000000]

setdiff(a,b)




# Beispiel-Arrays
array1 = [[i, i+1, i+2] for i in 1:10000000]
array2 = [[i+1000,i+1001,i+1002] for i in 1:100000]

# HAUPTPROBLEM?: [1, 3, 5] enspricht [5, 1, 3]

# Konvertiere Arrays in Sets für schnelles Nachschlagen
set1 = Set(array1)
set2 = Set(array2)

# Berechne die Differenzmenge
result_set = setdiff(set1, set2)

# Ausgabe des Ergebnisses
result = collect(result_set)


## #########

md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
print("tehrahedrons: ", length(md.Ω.faces))
print("swgfaces: ", length(md.X.fns))

@time swg_faces = swgfaces(md.Ω, md.Γ_nc, fast = true)
#@time swg_faces2 = swgfaces(md.Ω, md.Γ_nc, fast = false)
@time swg_faces_sa = ImpedancePredictionVIE.swgfaces_set_approach(md.Ω, md.Γ_nc)

##
using Test
# Andere Reihenfolge beim MultiTheading möglich mache daher normtest
@show sum(norm.(swg_faces))
@show sum(norm.(swg_faces2))
@test sum(norm.(swg_faces)) ≈ sum(norm.(swg_faces2))
@test sum(norm.(swg_faces)) ≈ sum(norm.(swg_faces3))