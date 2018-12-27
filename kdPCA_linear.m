function [X_d,X_star_d,X_hat_d,X_star_hat_d] = kdPCA_linear(X,X_gamma,lambda,X_star,R)

kappa = @(x,y) x*y';

if(nargout > 2)
    [X_d,X_star_d,X_hat_d,X_star_hat_d] = kdPCA(X,X_gamma,lambda,kappa,X_star,R);
else
    [X_d,X_star_d] = kdPCA(X,X_gamma,lambda,kappa,X_star,R);
end