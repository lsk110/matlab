function dn = cal_dn(Lap,r1,cn,dt)
cn_hat = fft2(cn);
dn_hat = cn_hat./(1 - dt*Lap.*(r1 + 1 + 2*Lap + Lap.^2));
dn = ifft2(dn_hat);

end