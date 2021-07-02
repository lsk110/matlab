function [rho_r] = aftertreat(psi,N,len)
rho_r = psi;
for i = 1:N
    for j = 1:N
        if i<len
        	ii = N-i;
            summ = sum(psi(ii:N,j)) ;
            summ1 =  sum(psi(1:i+len,j));
            summ = summ + summ1;
        end
        if i>N-len
            ii = i+len-N;
            summ = sum(psi(1:ii,j)) + sum(psi(i-len:N,j));
        end
        if i>=len && i<=N-len
            summ = sum(psi(i-len+1:i+len,j));
        end
        rho_r(i,j) = 1/len*summ;
    end
end