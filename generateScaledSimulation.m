function [X,X_gamma,X_star,T,L,L_star] = generateScaledSimulation(N)

T = 15;
C = 3;
C_star = 2;

M = T*C;
M2 = T*C_star;
% N = 100; 
sig  = 1;

W = randn(2,N);
L = zeros(M,2);
L_star = zeros(M2,2);

for ii = 1:C
    for jj = 1:T
        x = min(8,jj)*ii - 8/2*ii;
        y = max(0,jj-8)*ii - 8/2*ii;
        L(jj+(ii-1)*T,:) = [x y];
    end
end

for ii = 1:C_star
    ii2 = ii + 0.5;
    for jj = 1:T
        x = min(8,jj)*ii2 - 8/2*ii2;
        y = max(0,jj-8)*ii2 - 8/2*ii2;
        L_star(jj+(ii-1)*T,:) = [x y];
    end
end


X_0 = randn(M,N)*sig + (L*W);

[X,mu,st] = zscore(X_0);
X_star = (randn(M2,N)*sig + L_star*W - mu)./st;


X_gamma = zeros(M,N,3);

for ii = 1:C
    tts_ii = (1:T)+(ii-1)*T;
    X_gamma(tts_ii,:,2) = ones(length(tts_ii),1)*mean(X(tts_ii,:),1);
end
for ii = 1:T
    b1 = ii:T:M;
    X_gamma(b1,:,1) = ones(length(b1),1)*mean(X(b1,:),1);
end

X_gamma(:,:,3) = X-sum(X_gamma(:,:,1:2),3);
