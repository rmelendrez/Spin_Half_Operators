include("Sx_spin_half_Op.jl")
include("Sy_spin_half_Op.jl")
include("Sz_spin_half_Op.jl")

using LinearAlgebra
using .Sx_onsite_Mod
using .Sy_onsite_Mod
using .Sz_onsite_Mod




sites = 2
Na = sites
spin = 1/2
sps = 2
dim = Int64(2*spin+1)^Na
allsites = UnitRange(1,Na)
fullbasislist = UnitRange(0,(2^Na)-1)
fullbasislookup = Dict()
state = zeros(ComplexF64,dim)





# (↑ + i↓)(↑ + i↓)/2 = (↑↑ + i↑↓+ i↓↑ -↓↓)/2
state[1] = -1/2  # 00
state[2] = im/2 # 10
state[3] = im/2 # 01
state[4] = 1/2 # 11

expSx1 = Adjoint(state)*Sx_onsite_Mod.Sx_onsite(sites,1,spin)*state
expSy1 = Adjoint(state)*Sy_onsite_Mod.Sy_onsite(sites,1,spin)*state
expSz1 = Adjoint(state)*Sz_onsite_Mod.Sz_onsite(sites,1,spin)*state

println("Site 1 Spin expectation values")
println("<Sx1> = ", expSx1)
println("<Sy1> = ", expSy1)
println("<Sz1> = ", expSz1)


expSx2 = Adjoint(state)*Sx_onsite_Mod.Sx_onsite(sites,2,spin)*state
expSy2 = Adjoint(state)*Sy_onsite_Mod.Sy_onsite(sites,2,spin)*state
expSz2 = Adjoint(state)*Sz_onsite_Mod.Sz_onsite(sites,2,spin)*state

println("Site 2 Spin expectation value")
println("<Sx2> = ", expSx2)
println("<Sy2> = ", expSy2)
println("<Sz2> = ", expSz2)


expSx1Sx2 = Adjoint(state)*Sx_onsite_Mod.Sx_onsite(sites,1,spin)*Sx_onsite_Mod.Sx_onsite(sites,2,spin)*state
expSy1Sy2 = Adjoint(state)*Sy_onsite_Mod.Sy_onsite(sites,1,spin)*Sy_onsite_Mod.Sy_onsite(sites,2,spin)*state
expSz1Sz2 = Adjoint(state)*Sz_onsite_Mod.Sz_onsite(sites,1,spin)*Sz_onsite_Mod.Sz_onsite(sites,2,spin)*state



println("Two site correlation function")
println("<Sx1Sx2> = ", expSx1Sx2)
println("<Sy1Sy2> = ", expSy1Sy2)
println("<Sz1Sz2> = ", expSz1Sz2)







