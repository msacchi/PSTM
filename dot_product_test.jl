using PyPlot,SeisProcessing,SeisPlot,Printf

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



m1=randn(Nz,Nx)
d1=Operator_Kirk(m1, P, "forward")

d2=randn(Nt,Ntraces)
m2=Operator_Kirk(d2, P, "adj")

 dot_m = sum(m1.*m2) 
 dot_d = sum(d1.*d2)

 @printf "dot_m: %f\n" dot_m
 @printf "dot_d: %f\n" dot_d