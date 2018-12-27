function [X,X_gamma,X_star,T,L,L_star] = generateRotatedSimulation(N)

T = 9;
% angles = [0 pi pi*3/2 ];
% angles2 = pi/2;
angles = [0 -pi/2 pi  pi/2];
angles2 = [-pi/4 -pi*5/4];% pi*3/2];
% angles = [0 pi/2  ];
% angles2 = [pi];% pi*3/2];
C = length(angles);

M = T*C;
% N = 100; 
sig  = 1;

W = randn(2,N);
L = zeros(M,2);



% phis = linspace(pi,pi*3/2,T);

phis = linspace(deg2rad(10),deg2rad(80),T);

for ii = 1:C
    for jj = 1:T
        
        
        a = 1;
        b = 1;
        x = a*cos(phis(jj))*10;
        y = b*sin(phis(jj))*10;
%         
        x = max(0,8-jj) +1;
        y = max(0,jj-8) + 1;
        x = jj+1;
        y = 0;
        
        if(ii > 1)
            aa = angles(ii);
            R = [cos(aa) -sin(aa); sin(aa) cos(aa)];
            L(jj+(ii-1)*T,:) = L(jj,:)*R ;
        else
            L(jj+(ii-1)*T,:) = [x y]*0.25;
        end
    end
end
C2     = length(angles2);
L_star = zeros(C2*T,2);
for ii = 1:C2
    for jj = 1:T
        aa = angles2(ii);
        R = [cos(aa) -sin(aa); sin(aa) cos(aa)];
        L_star(jj+(ii-1)*T,:) = L(jj,:)*R ;
    end
end
M2 = size(L_star,1);



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