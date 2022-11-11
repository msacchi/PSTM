mutable struct Initialization_K

    # All the stuff I may need and julia programmers like to put
    # here to make things more complicated..

    dx::Float64
    dz::Float64
    dt::Float64
    Nx::Int64
    Nz::Int64
    Nt::Int64
    gx::Array{Float64,1}
    sx::Array{Float64,1}
    Ntraces::Int64
    x::Array{Float64,1}
    z::Array{Float64,1}
    c0::Float64
   
end


"""
    Simplified Forward and Adjoint Operators for Kirchhoff migration.

    Warning: No weights are applied, no proper aperture is designed,
    no interpolation is applied when converting traveltime to sample
    to grab from data.

    v is assumed constant but could be replaced by rms velocities and place
    the output at t0 (migration time) rather that z as I did below.

    - model,  m, is the reflectivity as a matrix m(1:Nz,1:Nx).

    - data are in format d(time,trace_number) and source and receiver
    positions are give as sx(trace_number) and gx(trace_number). So both
    input and output are matrices. time=1:Nt and trace_number=1:Ntraces.
    Ntrace = Ng*Ns (number of geophones x number of source), so I am assuming
    a fixed array of geophones that all record each moving source. A professional 
    version should get source and receiver positions from headers loop over total number
    of traces which probaly is not Ng*NS.

    - Still to do is also to divide ouput in offset bins.


    M D Sacchi 2022 after discussions with AQ and AA.

"""

function Operator_Kirk(In::Array{Float64,2}, P, Flag)
    
    if Flag=="adj"
        m = zeros(Float64,P.Nz,P.Nx)
        d = In
    else
        d=zeros(Float64,P.Nt,P.Ntraces)
        m = In
    end

    c0 = P.c0
    dt = P.dt
    
    if Flag=="adj"

        for  k = 1:P.Ntraces # Get trace
        sx = P.sx[k]
        gx = P.gx[k]
         for ix=1: P.Nx
            for iz=4: P.Nz
                d_s = sqrt((P.x[ix]-sx)^2+(P.z[iz])^2)
                d_g = sqrt((P.x[ix]-gx)^2+(P.z[iz])^2)
                t = (d_s+d_g)/c0
                it = convert(Int64,floor(t/dt))+1 # Interp is needed here
                if it<=P.Nt 
                    m[iz,ix] = m[iz,ix] + d[it,k] 
                end
            end
        end
    end

else

    for  k = 1:P.Ntraces # Get trace
        sx = P.sx[k]
        gx = P.gx[k]
         for ix =1: P.Nx
            for iz=4: P.Nz
                d_s = sqrt((P.x[ix]-sx)^2+(P.z[iz])^2)
                d_g = sqrt((P.x[ix]-gx)^2+(P.z[iz])^2)
                t = (d_s+d_g)/c0
                it = convert(Int64,floor(t/dt))+1 # Interp is needed here
                if it<=P.Nt 
                    d[it,k] = d[it,k] + m[iz,ix]
                end
            end
        end
    end

end

if Flag=="adj"
    Out = m
else
    Out = d
end

return Out
end