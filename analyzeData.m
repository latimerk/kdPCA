load('Data/kepecsProcessed.mat');

useFigureComposer = false;
saveResults = false;
figDir = 'tex/figs_src';

if(useFigureComposer)
    addpath ~/FigureComposer/matlab/
end

lengthScale = 50;
R = 20;

ms = 6;

fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;

if(~exist('DataFitting_RegularizationResults.mat','file'))
    lambdas   = linspace(0,5,1+100);
    lambdas = lambdas(2:end);
    % lambdas = [0 logspace(log10(1),log10(20e3),99)];
    L = length(lambdas);
    % lambdas_k = linspace(0,4,L);

    r2s_within_k = nan(L,1);
    r2s_within_d = nan(L,1);
    r2s_cv_k = nan(L,1);
    r2s_cv_d = nan(L,1);
    
    pv_d = nan(L,R,4);
    pv_k = nan(L,R,4);
    pv_star_d = nan(L,R,4);
    pv_star_k = nan(L,R,4);
    
    tse = sum(X(:).^2);
    tse_star = sum(X_star(:).^2);

    for ll = 1:L

        lambda_k = lambdas(ll);
        lambda   = lambdas(ll);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        C = size(X,1)/T;
        C_star = size(X_star,1)/T;

        if(ll == 1)
            results_kdPCA = struct();
            [results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d,K,K_star] = kdPCA_sqExp(X,X_gamma,lambda_k,lengthScale,X_star,R);
        else
            [results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA(X,X_gamma,lambda_k,K,X_star,R,K_star);
        end
        results_dPCA = struct();
        [ results_dPCA.X_d, results_dPCA.X_star_d, results_dPCA.X_hat_d, results_dPCA.X_star_hat_d] = dPCA(X,X_gamma,lambda,X_star,R);

        pv_star_k(ll,:,:) = 1-squeeze(sum(sum((X_star-results_kdPCA.X_star_hat_d).^2,1),2))./tse_star;
        pv_star_d(ll,:,:) = 1-squeeze(sum(sum((X_star-results_dPCA.X_star_hat_d ).^2,1),2))./tse_star;
        pv_k(ll,:,:)      = 1-squeeze(sum(sum((X     -results_kdPCA.X_hat_d     ).^2,1),2))./tse;
        pv_d(ll,:,:)      = 1-squeeze(sum(sum((X     -results_dPCA.X_hat_d      ).^2,1),2))./tse;
        

        %%
        rs = {results_kdPCA, results_dPCA};


        X_d = zeros(T,4,2,3,4);
        X_k = zeros(T,4,2,3,4);

        Xs_d = zeros(T,2,3,4);
        Xs_k = zeros(T,2,3,4);


        for dd = 1:2
            for cc = 1:4
                for rr = 1:3
                    for ii = 1:4
                        for jj = 1:2
                            if(dd == 1)
                                X_k(:,ii,jj,rr,cc) = rs{dd}.X_d(idxs{ii,jj},rr,cc);
                            elseif(dd == 2)
                                X_d(:,ii,jj,rr,cc) = rs{dd}.X_d(idxs{ii,jj},rr,cc);
                            end
                        end
                    end

                    if(dd == 1)
                        Xs_k(:,:,rr,cc) = reshape(rs{dd}.X_star_d(:,rr,cc),T,[]);
                    elseif(dd == 2)
                        Xs_d(:,:,rr,cc) = reshape(rs{dd}.X_star_d(:,rr,cc),T,[]);
                    end
                end
            end
        end





        ds_k = squeeze(mean(X_k(:,:,:,1,3),2));
        ds_d = squeeze(mean(X_d(:,:,:,1,3),2));


        r2s_d = zeros(5,2);
        r2s_k = zeros(5,2);

        r2s_d_avg = zeros(5,1);
        r2s_k_avg = zeros(5,1);

        for jj = 1:2
            for ii = 1:4
                ss = [1:(ii-1) (ii+1):4];
                ds_k_c = squeeze(mean(X_k(:,ss,:,1,3),2));
                ds_d_c = squeeze(mean(X_d(:,ss,:,1,3),2));

                r2s_d(ii,jj) = getR2(ds_d_c(:,jj),X_d(:,ii,jj,1,3));
                r2s_k(ii,jj) = getR2(ds_k_c(:,jj),X_k(:,ii,jj,1,3));
            end
            r2s_d(5,jj) = getR2(ds_d(:,jj),Xs_d(:,3-jj,1,3));
            r2s_k(5,jj) = getR2(ds_k(:,jj),Xs_k(:,3-jj,1,3));
        end
        for ii = 1:4
            ss = [1:(ii-1) (ii+1):4];
            ds_k_c = reshape(squeeze(mean(X_k(:,ss,:,1,3),2)),[],1);
            ds_d_c = reshape(squeeze(mean(X_d(:,ss,:,1,3),2)),[],1);

            r2s_d_avg(ii) = getR2(ds_d_c,reshape(squeeze(X_d(:,ii,:,1,3)),[],1));
            r2s_k_avg(ii) = getR2(ds_k_c,reshape(squeeze(X_k(:,ii,:,1,3)),[],1));
        end
        r2s_d_avg(5) = getR2(ds_d(:),reshape(Xs_d(:,end:-1:1,1,3),[],1));
        r2s_k_avg(5) = getR2(ds_k(:),reshape(Xs_k(:,end:-1:1,1,3),[],1));

        r2s_within_k(ll) = mean(r2s_k_avg(1:4));
        r2s_cv_k(ll)     = r2s_k_avg(5);
        r2s_within_d(ll) = mean(r2s_d_avg(1:4));
        r2s_cv_d(ll)     = r2s_d_avg(5);

        %%
        if(mod(ll,10) == 0 || ll == L)
            %%
            figure(10);
            clf
            subplot(1,3,1)
            hold on
            plot(lambdas,[r2s_within_k r2s_cv_k]);
            plot(lambdas,[r2s_within_d r2s_cv_d],'--');
            xlabel('regularization:','FontSize',fontSize_label);
            ylabel('R^2','FontSize',fontSize_label);
            set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
            hold off
            subplot(1,3,2)
            hold on
            plot(lambdas,[pv_k(:,1,3) pv_star_k(:,1,3)]);
            plot(lambdas,[pv_d(:,1,3) pv_star_d(:,1,3)],'--');
            xlabel('regularization:','FontSize',fontSize_label);
            ylabel('Variance explained (%)','FontSize',fontSize_label);
            title('First decision component','FontSize',fontSize_title);
            set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
            hold off
            subplot(1,3,3)
            hold on
            plot(lambdas,[sum(sum(pv_k(:,1:3,:),2),3) sum(sum(pv_star_k(:,1:3,:),2),3)]);
            plot(lambdas,[sum(sum(pv_d(:,1:3,:),2),3) sum(sum(pv_star_d(:,1:3,:),2),3)],'--');
            xlabel('regularization:','FontSize',fontSize_label);
            ylabel('Variance explained (%)','FontSize',fontSize_label);
            title('First three components (all factors)','FontSize',fontSize_title);
            set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
            hold off
            drawnow
        end
    end
    if(saveResults)
        save('DataFitting_RegularizationResults.mat','-v7.3','r2s_cv_k','r2s_cv_d','r2s_within_k','r2s_within_d','lambdas','lengthScale','pv_star_d','pv_star_k','pv_k','pv_d');
    end
    [~,ll_k] = max(r2s_cv_k);
    [~,ll_d] = max(r2s_cv_d);
    lambda = lambdas(ll_d);
    lambda_k = lambdas(ll_k);
    [results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d,K,K_star] = kdPCA(X,X_gamma,lambda_k,K,X_star,R,K_star);
    [ results_dPCA.X_d, results_dPCA.X_star_d, results_dPCA.X_hat_d, results_dPCA.X_star_hat_d] = dPCA(X,X_gamma,lambda,X_star,R);
else
    load('DataFitting_RegularizationResults.mat');
    [~,ll_k] = max(r2s_cv_k);
    [~,ll_d] = max(r2s_cv_d);
    lambda = lambdas(ll_d);
    lambda_k = lambdas(ll_k);
    
    [results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA_sqExp(X,X_gamma,lambda_k,lengthScale,X_star,R);
    [ results_dPCA.X_d, results_dPCA.X_star_d, results_dPCA.X_hat_d, results_dPCA.X_star_hat_d] = dPCA(X,X_gamma,lambda,X_star,R);
end
%%




%%
rs = {results_dPCA,results_kdPCA};

colors = [ 242 7 4;
          249 126 9;
          231 206 22;
          85  141 85;
          10 124 246;
          1  10 235]./255;
% colors = colors(end:-1:1,:);

X_d = zeros(T,4,2,3,4);
X_k = zeros(T,4,2,3,4);

Xs_d = zeros(T,2,3,4);
Xs_k = zeros(T,2,3,4);

NC = 2;
RR = length(rs);
figure(1);
clf;
tt_axis = (1:T)*0.01 - 0.5;
for dd = 1:RR
    
    for cc = 1:4
        for rr = 1:NC
            subplot(4,(NC+1)*RR-1,(cc-1)*((NC+1)*RR - 1) + rr + (NC+1)*(dd-1));
            hold on
            for ii = 1:4
                for jj = 1:2
                    if(jj == 1)
                        ls = '--';
                    else
                        ls = '-';
                    end
%                     plot(tt_axis,mean(rs{dd}.X_hat_d(idxs{ii,jj},:,rr,cc),2),ls,'Color',colors(ii+1,:));
                    plot(tt_axis,rs{dd}.X_d(idxs{ii,jj},rr,cc),ls,'Color',colors(ii+1,:))
                    if(dd == 1)
                        X_k(:,ii,jj,rr,cc) = rs{dd}.X_d(idxs{ii,jj},rr,cc);
                    elseif(dd == 2)
                        X_d(:,ii,jj,rr,cc) = rs{dd}.X_d(idxs{ii,jj},rr,cc);
                    end
                end


            end
%             plot(1:T,mean(rs{dd}.X_star_hat_d((1:T)+0*T,:,rr,cc),2),'-','Color',colors(1,:));
%             plot(1:T,mean(rs{dd}.X_star_hat_d((1:T)+1*T,:,rr,cc),2),'--' ,'Color',colors(6,:));
            plot(tt_axis,rs{dd}.X_star_d((1:T)+0*T,rr,cc),'-','Color',colors(1,:));
            plot(tt_axis,rs{dd}.X_star_d((1:T)+1*T,rr,cc),'--' ,'Color',colors(6,:));

            if(dd == 1)
                Xs_k(:,:,rr,cc) = reshape(rs{dd}.X_star_d(:,rr,cc),T,[]);
            elseif(dd == 2)
                Xs_d(:,:,rr,cc) = reshape(rs{dd}.X_star_d(:,rr,cc),T,[]);
            end
            hold off
            set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
            xlim([tt_axis(1) tt_axis(end)]);
        end
    end
end

drawnow;

set(gcf,'PaperSize',[15 12],'PaperPosition',[0 0 15 12]);

if(saveResults)
    if(useFigureComposer)
        matfig2fyp(gcf,sprintf('%s/Data_1_rough.fyp',figDir));
    else
        saveas(gcf,sprintf('%s/Data_1_rough.eps',figDir),'epsc');
    end
end


lw = 2;

cv_k_st    = '--';
cv_k_color = [0 0 0];
cv_d_st    = '--';
cv_d_color = [0.7 0.7 0.7];

tr_k_st    = '-';
tr_k_color = [0 0 0];
tr_d_st    = '-';
tr_d_color = [0.7 0.7 0.7];

figure(10);
clf

cc = 3;
rr = 1;
lw1 = 1;
lw2 = 2;
strs = {'dPCA','kdPCA'};
for dd = 1:2
    subplot(2,2,dd)
    hold on
    for jj = 1:2
            if(jj == 1)
                ls = '--';
            else
                ls = '-';
            end
        for ii = 1:4
            plot(tt_axis,rs{dd}.X_d(idxs{ii,jj},rr,cc),ls,'Color',colors(ii+1,:),'LineWidth',lw1)
        end
        if(dd == 1)
            plot(tt_axis,squeeze(mean(X_d(:,:,jj,rr,cc),2)),ls,'Color',[0 0 0],'LineWidth',lw2)
        elseif(dd == 2)
            plot(tt_axis,squeeze(mean(X_k(:,:,jj,rr,cc),2)),ls,'Color',[0 0 0],'LineWidth',lw2)
        end
        plot(tt_axis,rs{dd}.X_star_d((1:T)+0*T,rr,cc),'-','Color',colors(1,:),'LineWidth',lw2)
        plot(tt_axis,rs{dd}.X_star_d((1:T)+1*T,rr,cc),'--' ,'Color',colors(6,:),'LineWidth',lw2)

    end
    plot(tt_axis,rs{dd}.X_d(idxs{ii,jj},rr,cc),ls,'Color',colors(ii+1,:))
    set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
    title(sprintf('first decision component: %s',strs{dd}),'FontSize',fontSize_title);
    xlabel('time (s)','FontSize',fontSize_label);
    ylabel('projection','FontSize',fontSize_label);
    xlim([tt_axis(1) tt_axis(end)]);
    hold off
end

subplot(2,3,4)
hold on
plot(lambdas,r2s_within_k,tr_k_st,'LineWidth',lw,'Color',tr_k_color);
plot(lambdas,r2s_within_d,tr_d_st,'LineWidth',lw,'Color',tr_d_color);
plot(lambdas,r2s_cv_k,cv_k_st,'LineWidth',lw,'Color',cv_k_color);
plot(lambdas,r2s_cv_d,cv_d_st,'LineWidth',lw,'Color',cv_d_color);
set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
title('decision component stability','fontSize',fontSize_title);
xlabel('regularization:','FontSize',fontSize_label);
ylabel('R^2','FontSize',fontSize_label);
hold off
subplot(2,3,5)
hold on
plot(lambdas,pv_k(:,1,3)*100,tr_k_st,'LineWidth',lw,'Color',tr_k_color);
plot(lambdas,pv_d(:,1,3)*100,tr_d_st,'LineWidth',lw,'Color',tr_d_color);
plot(lambdas,pv_star_k(:,1,3)*100,cv_k_st,'LineWidth',lw,'Color',cv_k_color);
plot(lambdas,pv_star_d(:,1,3)*100,cv_d_st,'LineWidth',lw,'Color',cv_d_color);
set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
title('first decision component','FontSize',fontSize_title);
xlabel('regularization:','FontSize',fontSize_label);
ylabel('variance explained (%)','FontSize',fontSize_label);
hold off
subplot(2,3,6)
hold on
plot(lambdas,sum(sum(pv_k(:,1:NC,:),2),3)*100,tr_k_st,'LineWidth',lw,'Color',tr_k_color);
plot(lambdas,sum(sum(pv_d(:,1:NC,:),2),3)*100,tr_d_st,'LineWidth',lw,'Color',tr_d_color);
plot(lambdas,sum(sum(pv_star_k(:,1:NC,:),2),3)*100,cv_k_st,'LineWidth',lw,'Color',cv_k_color);
plot(lambdas,sum(sum(pv_star_d(:,1:NC,:),2),3)*100,cv_d_st,'LineWidth',lw,'Color',cv_d_color);
set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
title(sprintf('first %d components (all factors)',NC),'FontSize',fontSize_title);
xlabel('regularization:','FontSize',fontSize_label);
ylabel('variance explained (%)','FontSize',fontSize_label);
hold off

drawnow;

set(gcf,'PaperSize',[9 6],'PaperPosition',[0 0 9 6]);
if(saveResults)
    if(useFigureComposer)
        matfig2fyp(gcf,sprintf('%s/Data_2_rough.fyp',figDir));
    else
        saveas(gcf,sprintf('%s/Data_2_rough.eps',figDir),'epsc');
    end
end