%% 
load('Data/lipProcessed.mat');


H = 4;
baseline = nanmean(nanmean(nanmean(nanmean(X_1(:,Train_Cs,:,:,:),2),3),4),5);

X_star    = zeros( H*T,N);
X_star((1:T)+0*T,:) = squeeze(X_1(:,Test_Cs(1),1,1,:))';
X_star((1:T)+1*T,:) = squeeze(X_1(:,Test_Cs(2),2,1,:))';
X_star((1:T)+2*T,:) = squeeze(X_1(:,Test_Cs(1),1,2,:))';
X_star((1:T)+3*T,:) = squeeze(X_1(:,Test_Cs(2),2,2,:))';

% X_1  is (neuron, coherence, decision, targets, time)
% decision: 1=OUT,  2=IN
% targets: 1=false, 2=on



X_all  = X_1(:,Train_Cs,:,:,:);

X_star = X_star-baseline';
X_all = X_all-baseline;



C = size(X_all,2);
D = size(X_all,3);
E = size(X_all,4);
T = size(X_all,5);

X       = zeros( C*D*E*T,N);
X_gamma_0 = zeros( C*D*E*T,N,9);


X_t = mean(mean(mean(X_all,3),4),2);
X_s = mean(mean(mean(X_all,3),4),5);
X_d = mean(mean(mean(X_all,2),4),5);
X_f = mean(mean(mean(X_all,2),3),5);
X_dt = mean(mean(X_all - X_t - X_s - X_d - X_f,2),4);
X_st = mean(mean(X_all - X_t - X_s - X_d - X_f,3),4);
X_sd = mean(mean(X_all - X_t - X_s - X_d - X_f,4),5);
X_sf = mean(mean(X_all - X_t - X_s - X_d - X_f,3),5);
X_df = mean(mean(X_all - X_t - X_s - X_d - X_f,2),5);
X_ft = mean(mean(X_all - X_t - X_s - X_d - X_f,2),3);

X_dft = mean((X_all - X_t - X_s - X_d - X_f) - X_dt - X_ft  - X_df  - X_sf - X_sd - X_st,2);
X_sft = mean((X_all - X_t - X_s - X_d - X_f) - X_dt - X_ft  - X_df  - X_sf - X_sd - X_st,3);
X_sdt = mean((X_all - X_t - X_s - X_d - X_f) - X_dt - X_ft  - X_df  - X_sf - X_sd - X_st,4);

idxs = cell(C,D,E);

for ii = 1:C
    for jj = 1:D
        for kk = 1:E
            idxs{ii,jj,kk} =  (ii-1)*(D*T*E) + (jj-1)*(T*E) + (kk-1)*T + (1:T);
            for ll = 1:T
                rr = (ii-1)*(D*T*E) + (jj-1)*(T*E) + (kk-1)*T + ll;

                X(rr,:) = X_all(:,ii,jj,kk,ll)';

                X_gamma_0(rr,:,1) = X_t(:,1,1,1,ll);
                X_gamma_0(rr,:,2) = X_s(:,ii,1,1,1);
                X_gamma_0(rr,:,3) = X_d(:,1,jj,1,1);
                X_gamma_0(rr,:,4) = X_f(:,1,1,kk,1);

                X_gamma_0(rr,:,5) = X_st(:,ii,1,1,ll);
                X_gamma_0(rr,:,6) = X_dt(:,1,jj,1,ll);
                X_gamma_0(rr,:,7) = X_ft(:,1,1,kk,ll);

                X_gamma_0(rr,:,8) = X_df(:,1,jj,kk,1);
                X_gamma_0(rr,:,9) = X_dft(:,1,jj,kk,ll);
            end
        end
    end
end
X_gamma = zeros(C*D*T*E,N,5);
X_gamma(:,:,1) = X_gamma_0(:,:,1);
X_gamma(:,:,5) = sum(X_gamma_0(:,:,[2 5]),3);
X_gamma(:,:,4) = sum(X_gamma_0(:,:,[8 9]),3);
X_gamma(:,:,2) = sum(X_gamma_0(:,:,[3 6]),3);
X_gamma(:,:,3) = sum(X_gamma_0(:,:,[4 7]),3);

clear X_gamma_0;

save('Data/lipProcessed_marginals.mat','-v7.3','X_gamma','H','E','X','X_star','T','D','C','idxs','Test_Cs','Train_Cs');
