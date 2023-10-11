function geom()
 sx = zeros(20*30)
 gx = zeros(20*30)
 s0 = 10.0; smax = 980
 g0 = 15.0; gmax = 987
  k = 1
   for is=1:20
    for ig=1:30
        sx[k]=s0+(is-1)*(smax-s0)/(30-1)
        gx[k]=g0+(ig-1)*(gmax-g0)/(20-1)
        k=k+1
    end
  end
return sx,gx
end
