using BEAST
using CompScienceMeshes
using LinearAlgebra
using StaticArrays

##ihplates symmetric
src = CompScienceMeshes.meshsphere(1.0, 0.1)
#src2 = CompScienceMeshes.rotate(src, SVector(pi, 0, 0))

trg = CompScienceMeshes.translate(src, SVector(0, 1, -1))
Γsrc = src
Γtrg = trg


##
c = 3e8
f = 1e8
λ = c/f
k = 2*π/λ

MD = Maxwell3D.doublelayer(wavenumber=k)
MS = Maxwell3D.singlelayer(wavenumber=k)

Xsrc = raviartthomas(Γsrc)
Ytrg = buffachristiansen(Γtrg)

A = assemble(MD, Ytrg, Xsrc)

@views farblkasm = BEAST.blockassembler(
MD,
Ytrg,
Xsrc,
quadstrat=BEAST.defaultquadstrat(MD, Ytrg, Xsrc)
)

@views function farassembler(Z, tdata, sdata)
@views store(v,m,n) = (Z[m,n] += v)
farblkasm(tdata,sdata,store)
end

Z=zeros(ComplexF64, 10, 10)
farassembler(Z, Vector(1:10), Vector(1:10))
Z