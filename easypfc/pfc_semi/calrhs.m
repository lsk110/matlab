function [rhs] = calrhs(psi,dt)

    rhs = psi + dt*psi.^3;
    rhs = fft2(rhs);
end