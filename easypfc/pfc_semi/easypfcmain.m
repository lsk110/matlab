clear all;clc;close all;

N = 2^7;dt = 1e-10;

r = -0.45;
barpsi = -0.2;

[Dx,Dy,Lap,lap1,lap2] = operators(N);

% [x,y] = grid(N);

[LA] = coperator(Dx,Dy,Lap,lap1,lap2,dt,r);

[psi] = initialvalue(r,barpsi,N);

m = min(min(psi));
M = max(max(psi));

time = 0;

figure(1)
pcolor(psi);
colorbar
caxis([m,M]);
shading faceted
getframe;

energy = [];

for Nt = 1:10^3
    time = time + dt;
    if(time > 10)
        break;
    end
    [rhs] = calrhs(psi,dt);
    psi_hat = rhs./LA;
    psi = real(ifft2(psi_hat)) + dt*0.03*rand(N,N);
    
    if mod(Nt,1) == 0
        figure(1); 
		pcolor(psi);
        colorbar
        caxis([m,M]);
        shading interp;
        title(['T=',num2str(time)]);
        getframe;
    end
    
    energy1 = innerproduct(1/2*psi.*ifft2((r + 1 + 2*Lap + Lap.^2).*fft2(psi)) + psi.^4/4,N);
    energy = [energy,energy1];
end

figure(2)
plot(1e-4*(1:10^3),energy)