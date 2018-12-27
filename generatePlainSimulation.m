function [X,X_gamma,X_star,T,L,L_star] = generatePlainSimulation(N)

T = 15;
C = 3;

M = T*C;
% N = 100; 
sig  = 1;

W = randn(2,N);
L = zeros(M,2);

stimVector = [1 0.1]*5;
timeVector = [-1 0.2];

timePattern = (-7:7)'*timeVector;
stimPattern = (-1:1)'*stimVector;


for ii = 1:C
    L((1:T)+(ii-1)*T,:) = timePattern+stimPattern(ii,:);
end

L_star = 0.5*L((1:T)+(0)*T,:) + 0.5*L((1:T)+(1)*T,:);
L_star = [L_star;0.5*L((1:T)+(C-2)*T,:) + 0.5*L((1:T)+(C-1)*T,:)];

N2 = size(L_star,1);

X_0 = randn(M,N)*sig + (L*W);

[X,mu,st] = zscore(X_0);
X_star = (randn(N2,N)*sig + L_star*W - mu)./st;


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
