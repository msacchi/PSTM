using PyPlot,SeisProcessing,SeisPlot

include("pstm.jl")

    dx=10.0
    dz=10.0
    dt=0.004
    Nx=200
    Nz=240
    Nt=900

    c0 = 1500.0

    Ng = 100
    ## PLAY with number of sources to see improvement in migration
    Ns = 5

    x = zeros(Float64, Nx)
    z = zeros(Float64, Nz)

    for ix = 1:Nx
        x[ix]=(ix-1)*dx
    end
    for iz = 1:Nz
        z[iz]=(iz-1)*dz
    end
    
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

m[20:20+nw-1,20]=3*w
m[60:60+nw-1,140]=2*w
m[190:190+nw-1,150]=-2*w

m[50:50+nw-1,1:50].=-w
m[100:100+nw-1,51:end].=-w
m[70:70+nw-1,1:50].=-w
m[120:120+nw-1,51:end].=-w


# Forward d = L m

   d = Operator_Kirk(m, P, "foward")

# Adjoint

    ma= Operator_Kirk(d, P, "adj")

ma = ma/maximum(ma)
m = m/maximum(m)
d = d/maximum(d)
    figure(1)
    
    subplot(221); imshow(m,vmin=-0.3,vmax=0.3, cmap="seismic")
    title("True Image m")

    subplot(222); imshow(d,vmin=-0.3,vmax=0.3, cmap="seismic")
    title("Shots d = L . m")

    subplot(223); imshow(ma,vmin=-0.3,vmax=0.3, cmap="seismic")
    title("Migration ma = L' . d ")