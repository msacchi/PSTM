function Operator_Conv(In,adj; wavelet=ones(10))

   nt,nx = size(In)
   nw = length(wavelet)
   W = fft(vcat(wavelet,zeros(nt-nw))) 
   Out = zeros(nt,nx)

if (adj) 
  for k = 1:nx
    aux = In[:,k];  
    Aux = fft(aux).*conj.(W); 
    Out[:,k]  = real.(ifft(Aux))
  end
else
  for k = 1:nx
    aux = In[:,k];  
    Aux = fft(aux).*(W); 
    Out[:,k]  = real.(ifft(Aux))
  end
  end

return Out
end







