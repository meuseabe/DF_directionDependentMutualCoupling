function [c] = generate_MCM(M2)
cj = sqrt(-1);

for m = 1:M2
v(1,m) = ( (m^(-1) + m^(-1)*0.05.*randn(1,1)) + ...
        cj*(m^(-1) + m^(-1)*0.05.*randn(1,1)) )...
   /sqrt(2);
% v(1,m) = ( rand(1,1) + cj*(rand(1,1)) )/sqrt(2);
end

[vv,ii] = sort(abs(v),'descend'); c = transpose([1 v(ii(2:end))/2]);