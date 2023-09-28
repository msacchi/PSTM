using PyPlot, SeisProcessing, SeisPlot

include("pstm.jl")
#--------------

    dx = 10.0    
    dz = 10.0
    dt = 0.004
    Nx = 200
    Nz = 240
    Nt = 900
    c0 = 1500.0
    Ng = 100
    Ns = 5

    x = zeros(Float64, Nx)
    z = zeros(Float64, Nz)

    x = [(i-1)*dx for i=1:Nx]
    z = [(i-1)*dz for i=1:Nz]
  
    ds = 1.0*(Nx-1)*dx/Ns
    dg = 1.0*(Nx-1)*dx/Ng

    sx = zeros(Float64,Ng*Ns)
    gx = zeros(Float64,Ng*Ns)
  
    for is = 1: Ns
        for ig = 1: Ng
            j = (is-1)*Ng+ig
            sx[j] = (is-1)*ds
            gx[j] = (ig-1)*dg
        end
    end
    
    Ntraces = Ng*Ns
        
    P = Initialization_K(dx,dz,dt,Nx,Nz,Nt,gx,sx,Ntraces,x,z,c0)

m = zeros(Float64,Nz,Nx)
w = Ricker(dt=dz,f0=0.5/(2.0*dz)) # Wavelet in depth
nw = length(w)

# Simple earth model consisting of difractors and a vertical fault

m[50:50+nw-1,1:50].=-w
m[100:100+nw-1,51:end].=-w
m[70:70+nw-1,1:50].=-w
m[120:120+nw-1,51:end].=-w

  d = Operator_Kirk(m, P, "foward")
 ma = Operator_Kirk(d, P, "adj")

 mi = zeros(Float64,Nz,Nx)
 mu = 0.00001


