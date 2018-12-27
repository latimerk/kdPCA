%Calls kdPCA with the squared exponential (Gaussian) kernel with the given
%length scale.
%Arguments follow from kdPCA
function [X_d,X_star_d,X_hat_d,X_star_hat_d,K,K_star] = kdPCA_sqExp(X,X_gamma,lambda,lengthScale,X_star,R)

kappa = @(x,y) exp(-1/(2*lengthScale^2)*sum((x-y).^2,2));

if(nargout > 2)
    [X_d,X_star_d,X_hat_d,X_star_hat_d,K,K_star] = kdPCA(X,X_gamma,lambda,kappa,X_star,R);
else
    [X_d,X_star_d] = kdPCA(X,X_gamma,lambda,kappa,X_star,R);
end
    