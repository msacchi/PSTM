using PyPlot, SeisReconstruction, SeisProcessing, SeisPlot, FFTW

include("Operator_PSTM.jl")
include("Operator_Conv.jl")
include("geom.jl")

# Dot test product for 
# d  = W L m      Forward
# m' = L'W' d     Adjoint
# 
# where L: PST-Demigration Operator
#       W: Wavelet that is convolved with output of L


 wavelet = Ricker(); 

 Nx = 100
 Nz = 130
 dt = 0.004
 dx = 5.0
 dz = 5.0
 Nt = 400


 sx,gx = geom()
 Ntraces = length(sx) 

 Param_PSTM = Dict(:Ntraces=>Ntraces,:Nx=>Nx,:Nz=>Nz,:dt=>dt,
                   :dx=>dx,:dz=>dz,
                   :Nt=>Nt,:c0=>2000,:sx=>sx,:gx=>gx)

 Param_Conv = Dict(:wavelet=>wavelet)


 m1  = randn(Nz,Nx)
 d1 = Operator_Conv(Operator_PSTM(m1,false; Param_PSTM...), false; Param_Conv...)

 d2  = randn(Nt,Ntraces)
 m2 = Operator_PSTM(Operator_Conv(d2,true;  Param_Conv...), true;  Param_PSTM...)

 dot_d = sum(d1.*d2) 
 dot_m = sum(m1.*m2)

 r = abs(dot_d-dot_m)/(dot_d+dot_m)

 println( " ---------------------------------------------")
 println( " Passed the dot product test : ",   r< 1.0E-10)
 println( " ---------------------------------------------")
