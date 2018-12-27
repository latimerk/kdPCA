% dPCA by Kenneth Latimer (December 2018)
%
% Performs dPCA as described in
%  Kobak, Dmitry, et al. "Demixed principal component analysis of neural population data." Elife (2016)
%
%
% X is M x N matrix with N the number of Neurons, M is the total number of
%    observation
%
% X_gamma is a M x N x D tensor. Each X_gamma(:,:,i) is a matrix representing one experimental parameter for demixing.
%           This function assumes the use has already built each matrix and combined across any conditions.
%           sum(X_gamma,3) should be equal to X
%
% lambda is the regularization parameter
%
% R is the number of components for each parameter to compute
%
% X_star is a K x N matrix for test observations (not used for training)
%
% R is the number of components for each condition to compute
%
% X_d is a tensor  M x R x D,  X_g(i,j,d) is the X(i,:) projected onto the
%    jth component of parameter d
%
% X_hat_d is a tensor of M x N x R x D, where X_hat_d(i,:,j,d) is the
%    reconstroction of X(i,:) from the jth component of parameter d
%
% X_hat and X_star_hat_d are similar for the test conditions

function [X_d,X_star_d,X_hat_d,X_star_hat_d] = dPCA(X,X_gamma,lambda,X_star,R)

D = size(X_gamma,3);


N      = size(X,2);
M      = size(X,1);
M_star = size(X_star,1);


X_d      = nan(M,R,D);
X_star_d = nan(M_star,R,D);

if(nargout > 2)
    X_hat_d      = zeros(M,N,R,D);
    X_star_hat_d = zeros(M_star,N,R,D);
else
    X_hat_d = [];
    X_star_hat_d = [];
end


if(lambda > 0)
    lambda = lambda*norm(X,'fro')^2/M;
    A_inv = (X'*X + eye(N)*lambda)\X';
else
    A_inv = pinv(X'*X)*X';
end

for cc = 1:D
    
    A_k_ols = A_inv*X_gamma(:,:,cc);
    [uq_k_s,~,~] = svd(cov((X*A_k_ols)));
    H = uq_k_s(:,1:R);
    uk_s = A_k_ols*H;
    X_d(:,:,cc) = X*uk_s;
    
    X_star_d(:,:,cc) = X_star*uk_s;
    
    if(nargout > 2)
        dd = X*uk_s;
        dd_s = X_star*uk_s;

        for rr = 1:R
            X_hat_d(:,:,rr,cc) = dd(:,rr)*H(:,rr)';
            X_star_hat_d(:,:,rr,cc) = dd_s(:,rr)*H(:,rr)';
        end
    end
end

% 
% X_hat = sum(sum(X_hat_d(:,:,:,:),4),3);
% X_star_hat = sum(sum(X_star_hat_d(:,:,:,:),4),3);