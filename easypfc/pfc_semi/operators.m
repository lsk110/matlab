function [Dx,Dy,Lap,lap1,lap2] = operators(N)

k1 = [0:N/2-1 -N/2 -N/2+1:-1];
k2 = [0:N/2-1  0   -N/2+1:-1];

a1 = zeros(N,N); a2 = a1;
b1 = zeros(N,N); b2 = b1;

for i = 1:N
	a1(i,:) = k1 ;
	b1(:,i) = k1';
	a2(i,:) = k2 ;
	b2(:,i) = k2';
end

Dx = sqrt(-1) * a2;
Dy = sqrt(-1) * b2;

Lap = -(a1.^2 + b1.^2);


lap1 = -a1.^2+b1.^2;
lap2 = -a1.*b1;
end