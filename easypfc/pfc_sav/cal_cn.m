function cn = cal_cn(psi,dt,rn,bn,Lap,N)

bn_hat = fft2(bn);
bn_Lap = ifft2(Lap.*bn_hat);
cn = psi + dt*rn*bn_Lap - dt/2*innerproduct(bn.*psi,N).*bn_Lap;

end
