function [psi,x,y] = initialvalue(r,barpsi,N)
x = zeros(N*N,1);
y = zeros(N*N,1);
%triangular
a = 1;
At = 4/5*(barpsi + 1/3*sqrt(-15*r - 36*barpsi^2));
qt = 2*pi/a;
qt = sqrt(3)/2;

psi = zeros(N,N);

% for i = 1:N
%     for j = 1:N
for i = N/2-1:N/2+1
    for j = N/2-1:N/2+1
        psi(i,j) = At*(cos(qt*i)*cos(qt*j/sqrt(3)) - cos(2*qt*j/sqrt(3))/2) + barpsi;
%           psi(i,j) = 1;

%         if (i-N/2)^2 + (j-N/2)^2 < (N/5)^2
%             psi(i,j) = -0.21;
%         end
%         if abs(i-1/2*N)<1/4*N
%             psi(i,j) = barpsi;
%         end
    end
end

%%striped
% q = 1;
% A = sqrt(-4*(r/3 + barpsi^2));
% 
% psi = zeros(N,N);
% 
% for i = 1:N
%     for j = 1:N
%         x((i-1)*N+j) = i;
%         y((i-1)*N+j) = j;
%         psi (i,j) = A*sin(q*i) + barpsi;
%     end
% end



end