load DataFitting_RegularizationResults.mat
[~,ll_k] = max(r2s_cv_k);
[~,ll_d] = max(r2s_cv_d);

R = 2;


pv_d = squeeze(pv_d(ll_d,:,:));
pv_k = squeeze(pv_k(ll_k,:,:));

[~,rr_d] = sort(pv_d(:),'descend');
[~,rr_k] = sort(pv_k(:),'descend');



% rr_d = reshape(rr_d,[],4);
% rr_k = reshape(rr_k,[],4);

ss = {'independent','stim','decision','iteraction'};
for cc = 1:4
    fprintf('dPCA = %s\n',ss{cc});
    for ii = 1:R
        fprintf('\t Comp %d\n',ii);
        fprintf('\t    Var Exc: %.1f\n',pv_d(ii,cc)*100);
        
        ee = size(pv_d,1)*(cc-1)+ii;
        nn = find(rr_d == ee);
        
        fprintf('\t    Rank   : %d\n',nn);
    end
end

for cc = 1:4
    fprintf('kdPCA = %s\n',ss{cc});
    for ii = 1:R
        fprintf('\t Comp %d\n',ii);
        fprintf('\t    Var Exc: %.1f\n',pv_k(ii,cc)*100);
        
        ee = size(pv_k,1)*(cc-1)+ii;
        nn = find(rr_k == ee);
        
        fprintf('\t    Rank   : %d\n',nn);
    end
end