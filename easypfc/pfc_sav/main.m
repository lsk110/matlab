clear all;clc;close all;

N = 2^7;dt = 1e-9  ;

dt_min = 1e-9; dt_max = 1e-3; alphaalpha = 6e2;

r1 = -0.45;
barpsi = -0.2;
alpha = 0.7;

[Dx,Dy,Lap,lap1,lap2] = operators(N);

[psi,x,y] = initialvalue(r1,barpsi,N);

% psi = psi';
% data = [x,y,psi(:)];
% save('C:\Users\Administrator\Desktop\data.txt','data','-ascii');

m = min(min(psi));
M = max(max(psi));

time = 0;

figure(1)

subplot(1,2,1)
pcolor(psi);
% colorbar
caxis([m,M]);
shading interp
axis off
getframe;

% figure(2)
subplot(1,2,2)
rho_r = aftertreat(psi,N,10);
pcolor(rho_r);
caxis([m,M]);
shading interp
axis off
getframe;

energy = [];
ddt = [];

% filename1=strcat('psi',num2str(0));
% rear='.plt';
% filename1=strcat(filename1,rear);
% 
% fid1=fopen(filename1,'wt');
% % fprintf(fid,'VARIABLES=X,Y\n');
% % fprintf(fid, 'ZONE,N=%d,E=%d,F=FEPOINT,ET=TRIANGLE\n',N,N);
% fprintf(fid1, 'zone i=%d, j=%d, f=point\n',N,N);
% 
% % for i=1:1:N  
% %    
% %  fprintf(fid,'%-18.15f   %-18.15f    \n',VX(i),VY(i) );
% % end
% %   
% % 
% % for i=1:1:N
% %  fprintf(fid,'%d    %d    %d\n',EToV(i,1),EToV(i,2),EToV(i,3));
% % end
% for i = 1:N
%     for j = 1:N
%         fprintf(fid1,'%d    %d    %f\n',j,i,psi(i,j));
%     end
% end
%         
% fprintf(fid1,'\n');
% fclose(fid1);
% 
% filename2=strcat('rho',num2str(0));
% rear='.plt';
% filename2=strcat(filename2,rear);
% 
% fid2=fopen(filename2,'wt');
% % fprintf(fid,'VARIABLES=X,Y\n');
% % fprintf(fid, 'ZONE,N=%d,E=%d,F=FEPOINT,ET=TRIANGLE\n',N,N);
% fprintf(fid2, 'zone i=%d, j=%d, f=point\n',N,N);
% 
% % for i=1:1:N  
% %    
% %  fprintf(fid,'%-18.15f   %-18.15f    \n',VX(i),VY(i) );
% % end
%   
% 
% % for i=1:1:N
% %  fprintf(fid,'%d    %d    %d\n',EToV(i,1),EToV(i,2),EToV(i,3));
% % end
% for i = 1:N
%     for j = 1:N
%         fprintf(fid2,'%d    %d    %f\n',j,i,rho_r(i,j));
%     end
% end
%         
% fprintf(fid2,'\n');
% fclose(fid2);
% 
% 
% fid1=fopen(filename1,'a');
% fid2=fopen(filename2,'a');

for Nt = 1:10^5
    time = time + dt;
    if(time > 10)
        break;
    end
    
    rn = sqrt(innerproduct(psi.^4/4,N));
    bn = psi.^3/sqrt(innerproduct(psi.^4/4,N));
%     rn = sqrt(alpha)*rn;
%     bn = sqrt(alpha)*bn;
    cn = cal_cn(psi,dt,rn,bn,Lap,N);
    dn = cal_dn(Lap,r1,cn,dt);
    en = innerproduct(bn.*dn,N)/(1+dt*rn/2);
    psi = cal_dn(Lap,r1,1/2*dt*ifft2(Lap.*fft2(bn))*en + cn,dt);
    
    if mod(Nt,1) == 0
        figure(1); 
        subplot(1,2,1)
		pcolor(psi);
%         colorbar
        caxis([m,M]);
        axis off
        shading interp;
        title(['T=',num2str(time)]);
        getframe;
        
%         figure(2)
        subplot(1,2,2)
        rho_r = aftertreat(psi,N,10);
        pcolor(rho_r);
        caxis([m,M]);
        shading interp
        axis off
        getframe;
        
%         if mod(Nt,10) == 1
%         fprintf(fid1, 'zone i=%d, j=%d, f=point\n',N,N);
%         for i = 1:N
%             for j = 1:N
%                 fprintf(fid1,'%d    %d    %f\n',j,i,psi(i,j));
%             end
%         end
%         fprintf(fid1,'\n');
%         
%         fprintf(fid2, 'zone i=%d, j=%d, f=point\n',N,N);
%         for i = 1:N
%             for j = 1:N
%                 fprintf(fid2,'%d    %d    %f\n',j,i,rho_r(i,j));
%             end
%         end
%         fprintf(fid2,'\n');
%         end
    end
    energy1 = innerproduct(1/2*psi.*ifft2((r1 + 1 + 2*Lap + Lap.^2).*fft2(psi)) + psi.^4/4,N);
    energy = [energy,energy1];
    
    if Nt>1
        dt = max([dt_min,dt_max/sqrt((1+alphaalpha*abs(energy1-energy(end-1))/dt))]);
    end
    ddt = [ddt,dt];
end

fclose(fid1);
fclose(fid2);

figure(3)
b=1e-4*(1:10^3);
plot(b,energy)

figure(4)
plot(ddt);


