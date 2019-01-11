% kdPCA by Kenneth Latimer (December 2018)
%
% Performs kernel dPCA as described in
%  Latimer, Kenneth W. "Nonlinear demixed component analysis for neural population data as a low-rank kernel regression problem." arXiv preprint arXiv:1812.08238 (2018).
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
% kappa is the kernel function - k(x,y) must be able to accept x as a
% matrix where each row of the matrix is an observation and return a vector
% of the kernel between each row of x and the single vector in y
%
% NOTE: kappa can be replaced by a precomputed Gram matrix K and if K_star
%       is also precomputed. Otherwise, K_star is not needed.
%
% X_star is a K x N matrix for test observations (not used for training)
%
% R is the number of components for each parameter to compute
%
% X_d is a tensor  M x R x D,  X_g(i,j,d) is the X(i,:) projected onto the
%    jth component of parameter d
%
% X_hat_d is a tensor of M x N x R x D, where X_hat_d(i,:,j,d) is the
%    reconstroction of X(i,:) from the jth component of parameter d
%
% X_hat and X_star_hat_d are similar for the test conditions
%
% K is the Gram matrix (M x M)
% K_star is the K x M matrix of the kernel function between the training
%         and test data

function [X_d,X_star_d,X_hat_d,X_star_hat_d,K,K_star] = kdPCA(X,X_gamma,lambda,kappa,X_star,R,K_star)

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
    X_hat_d      = [];
    X_star_hat_d = [];
end

if(isa(kappa,'function_handle'))
    K = zeros(M,M);
    K_star = zeros(M_star,M);
    for ii = 1:M
        K(ii,ii:end) = kappa(X(ii,:),X(ii:end,:));
    end
    for ii = 1:M_star
        K_star(ii,:) = kappa(X_star(ii,:),X);
    end
    K = K+triu(K,1)';
else
    K = kappa;
end

if(lambda > 0)
    
    lambda = lambda*trace(K)./M;
    
    A_inv = (K + eye(M)*lambda);
else
    A_inv = pinv(K);
end

for cc = 1:D
    
    if(lambda > 0)
        A_k_ols = A_inv\X_gamma(:,:,cc);
    else
        A_k_ols = A_inv*X_gamma(:,:,cc);
    end
    [uq_k_s,~,~] = svd(cov((K*A_k_ols)));
    H = uq_k_s(:,1:R);
    uk_s = A_k_ols*H;

    
    X_d(:,:,cc) = K*uk_s;
    
    X_star_d(:,:,cc) = K_star*uk_s;
    
    if(nargout > 2)
        dd = K*uk_s;
        dd_s = K_star*uk_s;

        for rr = 1:R
            X_hat_d(:,:,rr,cc) = dd(:,rr)*H(:,rr)';
            X_star_hat_d(:,:,rr,cc) = dd_s(:,rr)*H(:,rr)';
        end
    end
end


% X_hat = sum(sum(X_hat_d(:,:,:,:),4),3);
% X_star_hat = sum(sum(X_star_hat_d(:,:,:,:),4),3);