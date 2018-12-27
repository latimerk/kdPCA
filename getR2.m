function [r2,r] = getR2(model,data)

r2 = 1 - sum((model-data).^2)./sum((data-mean(data)).^2);

[aa,~] = corrcoef(model,data);
r = aa(1,2);