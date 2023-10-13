"""
Operator_Conv: Convolution of traces with a wavelet when adj=false
               Crosscorrelation of traces a wavelet when adj=true 
"""

function Operator_Conv(In,adj; wavelet=ones(10))

   nt,nx = size(In)
   nf = nextpow(2,nt)
   nw = length(wavelet)
   W = fft(vcat(wavelet,zeros(nf-nw))) 
   Out = zeros(nt,nx)

if (adj) 
  for k = 1:nx
    aux = In[:,k];  
    Aux = fft(vcat(aux,zeros(nf-nt))).*conj.(W); 
    tmp = real.(ifft(Aux))
    Out[:,k]  = tmp[1:nt] 
  end
else
  for k = 1:nx
    aux = In[:,k];  
    Aux = fft(vcat(aux,zeros(nf-nt))).*(W); 
    tmp = real.(ifft(Aux))
    Out[:,k]  = tmp[1:Nt] 
  end
  end

return Out
end







