using PyPlot, SeisReconstruction, SeisProcessing, SeisPlot, FFTW
include("Operator_PSTM.jl")
include("Operator_Conv.jl")
include("geom.jl")

# Inverse problem crime. Forward demigration operator is used to generate data and
# and then CGLS is used to solve the LSM (Least-squares migration) problem

 wavelet = Ricker(dt=0.004, f0=15) 

 Nx = 800
 Nz = 530
 dt = 0.004
 dx = 5.0
 dz = 5.0
 Nt = 600

 sx,gx = geom()

 Ntraces = length(sx) 

 Param_PSTM = Dict(:Ntraces=>Ntraces,:Nx=>Nx,:Nz=>Nz,:dt=>dt,
                   :dx=>dx,:dz=>dz,
                   :Nt=>Nt,:c0=>3000,
                   :sx=>sx,:gx=>gx)

 Param_Conv = Dict(:wavelet=>wavelet)


 m = zeros(Nz,Nx);
 m[140,140] =1;
 m[200,260:end] .=-1;
 m[160,40:600] .=-1;
 m[160,560] =1;

# Data

 d  = Operator_Conv(Operator_PSTM(m,false; Param_PSTM...), false; Param_Conv...)

 # Migration 
mig = Operator_PSTM(Operator_Conv(d,true; Param_Conv...), true; Param_PSTM...)


   # LS Migration.. see how I pass the operators 
  
x,J = ConjugateGradients(d,[Operator_Conv, Operator_PSTM],[Param_Conv, Param_PSTM]; Niter=40,mu=0.01)


dpred  = Operator_Conv(Operator_PSTM(x,false; Param_PSTM...), false; Param_Conv...)

SeisPlotTX(d/maximum(abs.(d)), vmin=-1, vmax=1, cmap="seismic", fignum=1, hbox=4, wbox=8, 
                               dx=1, dy=dt, xlabel="Trace number", ylabel="t (s)", title="Data")
tight_layout()
savefig("figure1.png")

SeisPlotTX(mig/maximum(abs.(mig)), vmin=-1, vmax=1, cmap="seismic", fignum=2, aspect="equal", hbox=4, wbox=8, 
                               dx=dx, dy=dz, xlabel="x (m)", ylabel="z (m)", title="PSTM")
tight_layout()
savefig("figure2.png")


SeisPlotTX(x/maximum(abs.(x)), vmin=-1, vmax=1, cmap="seismic", fignum=3, aspect="equal", hbox=4, wbox=8, 
                               dx=dx, dy=dz, xlabel="x (m)", ylabel="z (m)", title="LS-PSTM")
tight_layout()
savefig("figure3.png")


figure(4)
plot(J);ylabel("Normalized Cost"); xlabel("Iteration")

tight_layout()
savefig("figure4.png")
