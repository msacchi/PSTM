"""
Operator_PSTM: Forward and Adjoint operator for PSTM
               Constant velocity PSTM
"""

function Operator_PSTM(In, adj; Ntraces=20, Nx=10, Nz=10, dt=0.004, 
                                Nt=100, c0=2000, sx=rand(10), gx=rand(10),
                                dx=10, dz=10)

   
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
                tt = t/dt+1.0
                it1 = convert(Int64,floor(tt))
                alpha = tt  - it1*1.0;
                it2 = it1+1
                if it2<= Nt
                    m[iz,ix] = m[iz,ix] +  (1.0-alpha)*d[it1,k]+alpha*d[it2,k]
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
                tt = t/dt+1.0
                it1 = convert(Int64,floor(tt))
                alpha = tt  - it1*1.0;
                it2 = it1+1
                if it2<= Nt
                    d[it1,k] = d[it1,k] + (1.0-alpha)*m[iz,ix]
                    d[it2,k] = d[it2,k] +       alpha*m[iz,ix]
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


