function Operator_PSTM(In, adj; Ntraces=20, Nx=10, Nz=10, dt=0.004, 
                    Nt=100, c0=2000, sx=rand(10), gx=rand(10),
                    dx=10, dz=10)

# Simple constant v PSTM 
   
    if (adj)
        m = zeros(Float64,Nz,Nx)
        d = In
    else
        d=zeros(Float64,Nt,Ntraces)
        m = In
    end

 x = collect(0:dx:(Nx-1)*dx)
 z = collect(0:dz:(Nz-1)*dz)
   
    if (adj)
        
        for  k = 1:Ntraces # Get trace
            for ix = 1:Nx
                for iz = 4:Nz
                d_s = sqrt((x[ix]-sx[k])^2+z[iz]^2)
                d_g = sqrt((x[ix]-gx[k])^2+z[iz]^2)
                t = (d_s+d_g)/c0
                it = convert(Int64,floor(t/dt))+1 # Interp is needed here
                if it<= Nt
                    m[iz,ix] = m[iz,ix] + d[it,k]
                end
                end
            end
        end
        
    else

     for  k = 1:Ntraces # Get trace
            for ix = 1:Nx
                for iz = 4:Nz
                d_s = sqrt((x[ix]-sx[k])^2+z[iz]^2)
                d_g = sqrt((x[ix]-gx[k])^2+z[iz]^2)
                t = (d_s+d_g)/c0
                it = convert(Int64,floor(t/dt))+1 # Interp is needed here
                if it<= Nt
                    d[it,k] = d[it,k] + m[iz,ix]
                end
                end
            end
        end
    end
    
    
    
if (adj)
    Out = m
else
    Out = d
end

return Out
end


