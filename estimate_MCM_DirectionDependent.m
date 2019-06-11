function [c Call] = estimate_MCM_DirectionDependent(g,Ab,U,ZC)
ZC = ZC';
M2_C = size(ZC,1);
[M,K] = size(Ab);
G = [];
TT = [];
for k = 1:K
    for m = 1:M2_C
        Tc(:,1*(m-1)+1:1*m) = ZC(m,1).E*Ab(:,k);
    end
    TT = [TT Tc];
    %[Tc] = transform_MCM(M2_C,diag(g)*Ab(:,k),ZC);
    %P((M-K)*(k-1)+1:(M-K)*k,1:M2_C) = U'*Tc;
    
    GG(k,:,:) = Tc'*(U*U')*Tc;
    G = blkdiag(G,squeeze(GG(k,:,:)));
end
%% solve with cvx
% Gr = real(G);
% Gi = imag(G);
% cvx_begin quiet
% %cvx_solver SeDuMi
% variable cr(M2_C*K)
% variable ci(M2_C*K)
% minimize ( trace(([cr;ci]'*[Gr Gi;-Gi Gr]*[cr;ci])) )
% subject to
% 
% cr(1 + ([1:K]-1)*M2_C,1) == 1
% ci(1 + ([1:K]-1)*M2_C,1) == 0
% cvx_end
% 
% 
% c1 = cr - 1i*ci;
%%
% 
% % % Method of Friedlander.
% M2_C = M2_C;
W = zeros(M2_C,1); u = zeros(1,1);
W(1,1) = 1; u(1,1) = 1;
for k = 1:K
    c(1 + (k-1)*M2_C:k*M2_C,1) = inv(squeeze(GG(k,:,:)))*W*inv(W'*inv(squeeze(GG(k,:,:)))*W)*u;
end
% [c1 c]
% %%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:K
C = zeros(M,M);
for m = 1:M2_C;
    C = C +  ZC(m,1).E*c(M2_C*(k-1) + m,1);
end
Call(k,:,:) = C;
end



