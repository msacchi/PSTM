"""
geom: a simple function to set the geometry for a PSTM experiment
"""

function geom()

 sx = zeros(20*70)
 gx = zeros(20*70)
 s0 = 30.0; smax = 3800
 g0 = 45.0; gmax = 3500
  k = 1
   for is=1:20
    for ig=1:70
        sx[k]=s0+(is-1)*(smax-s0)/(30-1)
        gx[k]=g0+(ig-1)*(gmax-g0)/(70-1)
        k=k+1
    end
  end
return sx,gx
end
