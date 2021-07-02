function [psi] = initialvalue(r,barpsi,N)

% %triangular
% a = 1;
% At = 4/5*(barpsi + 1/3*sqrt(-15*r - 36*barpsi^2));
% qt = 2*pi/a;
% qt = sqrt(3)/2;
% 
% psi = zeros(N,N);
% 
% for i = 1:N
%     for j = 1:N
%         psi(i,j) = At*(cos(qt*i)*cos(qt*j/sqrt(3)) - cos(2*qt*j/sqrt(3))/2) + barpsi;
%     end
% end

%%striped
q = 1;
A = sqrt(-4*(r/3 + barpsi^2));

psi = zeros(N,N);

for i = 1:N
    for j = 1:N
        psi (i,j) = A*sin(q*i) + barpsi;
    end
end



end