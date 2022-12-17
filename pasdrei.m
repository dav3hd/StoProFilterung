function [p] = pasdrei(m)

M = m+1;
p = eye(M);
p(2:M,1) = 1.0;
p(3:M,2) = 2:M-1;

for n=4:M 
    for k=3:M-1
        p(n,k) = p(n-1,k-1)+p(n-1,k);
    end
end


end

