function [LA] = coperator(Dx,Dy,Lap,lap1,lap2,dt,r)
LA = 1 + dt*(-(r+1)*Lap - 2*Lap.^2 - Lap.^3);
end