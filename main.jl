using PyPlot, SeisReconstruction, SeisProcessing, SeisPlot, FFTW

include("Operator_PSTM.jl")
include("Operator_Conv.jl")
include("geom.jl")

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


m = zeros(Nz,Nx);
m[40,40] =1;
m[60,60] =-1;

d  = Operator_Conv(Operator_PSTM(m,false; Param_PSTM...), false; Param_Conv...)

Param_Conv = Dict(:wavelet=>wavelet)
Operator = [Operator_Conv, Operator_PSTM]
  Param  = [Param_Conv, Param_PSTM]

dpred  = Operator_Conv(Operator_PSTM(x,false; Param_PSTM...), false; Param_Conv...)
x,J = ConjugateGradients(d,Operator,Param; Niter=130,mu=0.01)
subplot(221);SeisPlotTX(x/maximum(abs.(x)),vmin=-1,vmax=1;cmap="seismic",fignum=1)
tight_layout()
