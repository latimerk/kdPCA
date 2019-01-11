function [ResultsSummary] = simulationFigures(generator,includeLinear,N_sims)

lambda_k = 1;
lambda_kl = 1;
lambda   = 1;
lengthScale = 5;


R = 2;

if(nargin < 3)
    N_sims = 100^2;
end

ms = 6;
lw  = 1.5;
lw2 = 1;

fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,X_gamma,X_star,T,L,L_star] = generator();
% MM = cov(X)+eye(size(X,2))*4;

C = size(X,1)/T;
C_star = size(X_star,1)/T;

results_kdPCA = struct();
[results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA_sqExp(X,X_gamma,lambda_k,lengthScale,X_star,R);
% [results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA_md(X,X_gamma,lambda_k,MM,X_star,R);

if(includeLinear)
    results_lkdPCA = struct();
    [results_lkdPCA.X_d,results_lkdPCA.X_star_d,results_lkdPCA.X_hat_d,results_lkdPCA.X_star_hat_d] = kdPCA_linear(X,X_gamma,lambda_kl,X_star,R);
end

results_dPCA = struct();
[ results_dPCA.X_d, results_dPCA.X_star_d, results_dPCA.X_hat_d, results_dPCA.X_star_hat_d] = dPCA(X,X_gamma,lambda,X_star,R);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nm = 2;
if(includeLinear)
    nm = 3;
    rs = {results_dPCA,results_lkdPCA,results_kdPCA};
else
    rs = {results_dPCA,results_kdPCA};
end
figure(1);
clf
subplot(nm+1,1,1);
hold on
colors_c = zeros(C,3);
colors_c_star = zeros(C_star,3);
for ii = 1:C
    tts = (1:T) + (ii-1)*T;
    pp = plot(L(tts,1),L(tts,2),'o-','LineWidth',lw,'MarkerSize',ms);
    colors_c(ii,:) = pp.Color;
end
for ii = 1:C_star
    tts = (1:T) + (ii-1)*T;
    pp = plot(L_star(tts,1),L_star(tts,2),'x--','LineWidth',lw,'MarkerSize',ms);
    colors_c_star(ii,:) = pp.Color;
end
for ii = 1:C
    tts = (1:T) + (ii-1)*T;
    plot(L(tts(1),1),L(tts(1),2),'o-','Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
end
for ii = 1:C_star
    tts = (1:T) + (ii-1)*T;
    plot(L_star(tts(1),1),L_star(tts(1),2),'o--','Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
end
title('true latent space','FontSize',fontSize_title);
xlabel('latent dim. 1','FontSize',fontSize_label);
ylabel('latent dim. 2','FontSize',fontSize_label);
set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
axis equal
axis square
hold off

Ts       = nan(T,C     ,nm);
Ts_star  = nan(T,C_star,nm);
mus      = nan(C       ,nm);
mus_star = nan(C_star  ,nm);
sig      = nan(C       ,nm);
sig_star = nan(C_star  ,nm);

for cc = 1:(nm)
    for jj = 1:3
        subplot(nm+1,3,(cc*3)+jj);
        hold on
        
        
        for ii = 1:C
            tts = (1:T) + (ii-1)*T;
            plot(rs{cc}.X_d(tts   ,1,jj),rs{cc}.X_d(tts   ,2,jj),'o-','Color',colors_c(ii,:),'LineWidth',lw,'MarkerSize',ms);
            plot(rs{cc}.X_d(tts(1),1,jj),rs{cc}.X_d(tts(1),2,jj),'o-','Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
            
            if(jj == 2)
                mus(ii,cc) = mean(rs{cc}.X_d(tts   ,1,jj));
                sig(ii,cc) = std(rs{cc}.X_d(tts   ,1,jj));
            end
        end
        for ii = 1:C_star
            tts = (1:T) + (ii-1)*T;
            plot(rs{cc}.X_star_d(tts   ,1,jj),rs{cc}.X_star_d(tts   ,2,jj),'x--','Color',colors_c_star(ii,:),'LineWidth',lw,'MarkerSize',ms);
            plot(rs{cc}.X_star_d(tts(1),1,jj),rs{cc}.X_star_d(tts(1),2,jj),'o--','Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
        
            if(jj == 2)
                mus_star(ii,cc) = mean(rs{cc}.X_star_d(tts   ,1,jj));
                sig_star(ii,cc) = std(rs{cc}.X_star_d(tts   ,1,jj));
            end
        end
        hold off
        
        if(cc == 1)
            if(jj == 1)
                title('estimated time space','FontSize',fontSize_title);
            elseif(jj == 2)
                title('estimated stimulus space','FontSize',fontSize_title);
            else
                title('estimated interaction space','FontSize',fontSize_title);
            end
        end
        xlabel('component 1','FontSize',fontSize_label);
        ylabel('component 2','FontSize',fontSize_label);
        axis equal
        set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    end
    
    Ts(:,:,cc)      = reshape(rs{cc}.X_d(:,1,jj),T,[]);
    Ts_star(:,:,cc) = reshape(rs{cc}.X_star_d(:,1,jj),T,[]);
end

set(gcf,'PaperSize',[9 (nm+1)*3],'PaperPosition',[0 0 9 (nm+1)*3]);
pause(1);
drawnow;
pause(1);


%%

figure(3);
clf



for cc = 1:(nm)
    for jj = 1:4
        subplot(nm+1,4,(cc*4)+jj);
        hold on
        
        switch jj
            case 1
                ccs = [1 3];
                dds = [1 2;
                       1 2];
            case 2
                ccs = [2 3];
                dds = [1 2;
                       1 2];
            case 3
                ccs = [1 2];
                dds = [1 2;
                       1 2];
            case 4
                ccs = [1 2 3];
                dds = [1 2;
                       1 2;
                       1 2];
        end
        
        for ii = 1:C
            tts = (1:T) + (ii-1)*T;
            xx = zeros(length(tts),2);
            
            for kk = 1:length(ccs)
                xx(:,1) = xx(:,1) + rs{cc}.X_d(tts   ,dds(kk,1),ccs(kk));
                xx(:,2) = xx(:,2) + rs{cc}.X_d(tts   ,dds(kk,2),ccs(kk));
            end
            
            plot(xx(:,1),xx(:,2),'o-','Color',colors_c(ii,:),'LineWidth',lw,'MarkerSize',ms);
            plot(xx(1,1),xx(1,2),'o-','Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
            
        end
        for ii = 1:C_star
            tts = (1:T) + (ii-1)*T;
            xx = zeros(length(tts),2);
            for kk = 1:length(ccs)
                xx(:,1) = xx(:,1) + rs{cc}.X_star_d(tts   ,dds(kk,1),ccs(kk));
                xx(:,2) = xx(:,2) + rs{cc}.X_star_d(tts   ,dds(kk,2),ccs(kk));
            end
            plot(xx(:,1),xx(:,2),'x--','Color',colors_c_star(ii,:),'LineWidth',lw,'MarkerSize',ms);
            plot(xx(1,1),xx(1,2),'o--','Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
        
        end
        hold off
        
        if(cc == 1)
            if(jj == 1)
                title('time + inter','FontSize',fontSize_title);
            elseif(jj == 2)
                title('stim + inter','FontSize',fontSize_title);
            elseif(jj == 3)
                title('time + stim','FontSize',fontSize_title);
            else
                title('time + stim + inter','FontSize',fontSize_title);
            end
        end
        xlabel('component 1','FontSize',fontSize_label);
        ylabel('component 2','FontSize',fontSize_label);
        axis equal
        set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    end
    
end
%%
set(gcf,'PaperSize',[12 (nm+1)*3],'PaperPosition',[0 0 12 (nm+1)*3]);
pause(1);
drawnow;
pause(1);
return;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dp_s = nan(N_sims,nm);
dp_s_star = nan(N_sims,nm);
CC_t = nan(N_sims,nm);
CC_t_star = nan(N_sims,nm);

tt_c      = repmat((1:T)',1,C);
tt_c_star = repmat((1:T)',1,C_star);

perVar = nan(N_sims,3,nm);
perVar_star = nan(N_sims,3,nm);

parfor ss = 1:N_sims
    if(mod(ss,50) == 0)
        fprintf('sim %d / %d\n',ss,N_sims);
    end
    [X,X_gamma,X_star] = generator();
    
    pv = zeros(3,nm);
    pvs = zeros(3,nm);
    
    tse = sum(X(:).^2);
    tse_star = sum(X_star(:).^2);

    results_kdPCA_c = struct();
    [results_kdPCA_c.X_d,results_kdPCA_c.X_star_d,results_kdPCA_c.X_hat_d,results_kdPCA_c.X_star_hat_d] = kdPCA_sqExp(X,X_gamma,lambda_k,lengthScale,X_star,R);

    results_lkdPCA_c = struct();
    if(includeLinear)
        [results_lkdPCA_c.X_d,results_lkdPCA_c.X_star_d,results_lkdPCA_c.X_hat_d,results_lkdPCA_c.X_star_hat_d] = kdPCA_linear(X,X_gamma,lambda_kl,X_star,R);
    end

    results_dPCA_c = struct();
    [ results_dPCA_c.X_d, results_dPCA_c.X_star_d, results_dPCA_c.X_hat_d, results_dPCA_c.X_star_hat_d] = dPCA(X,X_gamma,lambda,X_star,R);

    if(includeLinear)
        rs_c = {results_dPCA_c,results_lkdPCA_c,results_kdPCA_c};
    else
        rs_c = {results_dPCA_c,results_kdPCA_c};
    end
    
    for ii = 1:nm
        
        Ts_c      = reshape(rs_c{ii}.X_d(:,1,1),T,[]);
        Ts_star_c = reshape(rs_c{ii}.X_star_d(:,1,1),T,[]);
        
        [aa,~] = corrcoef(Ts_c(:),tt_c(:));
        CC_t(ss,ii) = aa(1,2);
        
        
        [aa,~] = corrcoef([Ts_star_c(:);Ts_c(:)],[tt_c_star(:);tt_c(:)]);
        CC_t_star(ss,ii) = aa(1,2);
        
        
        dps = nan(C,C+C_star);
        
        for jj = 1:C
            A = rs_c{ii}.X_d((1:T)+(jj-1)*T,1,2);
            for kk = (jj+1):C
                B = rs_c{ii}.X_d((1:T)+(kk-1)*T,1,2);
                dps(jj,kk) = dPrime(A,B);
            end
            for kk = 1:C_star
                B = rs_c{ii}.X_star_d((1:T)+(kk-1)*T,1,2);
                dps(jj,C+kk) = dPrime(A,B);
            end
        end
        dp_s(ss,ii) = nanmin(nanmin(abs(dps(1:C,1:C))));
        dp_s_star(ss,ii) = nanmin(nanmin(abs(dps(:,:))));
        
        for jj = 1:3
            se = (X-rs_c{ii}.X_hat_d(:,:,1,jj)).^2;
            sse = sum(se(:));
            pv(jj,ii) = 1-sse./tse;
            
            se = (X_star-rs_c{ii}.X_star_hat_d(:,:,1,jj)).^2;
            sse = sum(se(:));
            pvs(jj,ii) = 1-sse./tse_star;
        end
    end
    perVar(ss,:,:) = pv;
    perVar_star(ss,:,:) = pvs;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
clf

CC_t = CC_t.^2;
CC_t_star = CC_t_star.^2;


for cc = 1:nm
    %%
    subplot(nm,5,(cc-1)*5 + 1);
    hold on
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        plot(1:T,rs{cc}.X_d(tts   ,1,1),'o-','Color',colors_c(ii,:),'LineWidth',lw,'MarkerSize',ms);
        plot(1  ,rs{cc}.X_d(tts(1),1,1),'o-','Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        plot(1:T,rs{cc}.X_star_d(tts   ,1,1),'x--','Color',colors_c_star(ii,:),'LineWidth',lw,'MarkerSize',ms);
        plot(1  ,rs{cc}.X_star_d(tts(1),1,1),'o--','Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'LineWidth',lw2,'MarkerSize',ms);
    end
    ylabel('time PC 1','FontSize',fontSize_label);
    xlabel('time','FontSize',fontSize_label);
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    hold off
    
    
    %%
    
    ctr_min = min([abs(CC_t(:));abs(CC_t_star(:))]);
    ctr_max = max([abs(CC_t(:));abs(CC_t_star(:))]);
    
    ctr_min = floor(ctr_min*10)/10;
    ctr_max = ceil(ctr_max*10)/10;
    
    ctrs = linspace(ctr_min,ctr_max,31);
    h1 = hist(abs(CC_t(:,cc)),ctrs);
    h2 = hist(abs(CC_t_star(:,cc)),ctrs);
%     h1 = h1./sum(h1);
%     h2 = h2./sum(h2);
    
    subplot(nm,5,(cc-1)*5 + 2);
    hold on
    
    plot(ctrs,h1,'Color',[0 0 0],'LineWidth',lw);
    plot(ctrs,h2,'--','Color',[0.5 0.5 0.5],'LineWidth',lw);
    
    xlabel('correlation','FontSize',fontSize_label);
    ylabel('n. sims','FontSize',fontSize_label);
    xlim([ctrs(1) ctrs(end)]);
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    hold off
    
    
    
    %%
    subplot(nm,5,(cc-1)*5 + 3);
    hold on
    
    mu_min = min([mus(:,cc)-sig(:,cc)*3;mus_star(:,cc)-sig_star(:,cc)*3]);
    mu_max = max([mus(:,cc)+sig(:,cc)*3;mus_star(:,cc)+sig_star(:,cc)*3]);
    
    tts_dist = linspace(mu_min,mu_max,100);
    
    for ii = 1:C
        plot(tts_dist,normpdf(tts_dist,mus(ii,cc),sig(ii,cc)),'Color',colors_c(ii,:),'LineWidth',lw);
    end
    for ii = 1:C_star
        plot(tts_dist,normpdf(tts_dist,mus_star(ii,cc),sig_star(ii,cc)),'--','Color',colors_c_star(ii,:),'LineWidth',lw);
    end
    xlim([tts_dist(1), tts_dist(end)]);
    ylabel('approx. frequency','FontSize',fontSize_label);
    xlabel('stimulus PC 1','FontSize',fontSize_label);
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    hold off
    
    
    %%
    
%     ctrs = linspace(min([dp_s(:);dp_s_star(:)]),max([dp_s(:);dp_s_star(:)]), 31);
%     h1 = hist(dp_s(:,cc),ctrs);
%     h2 = hist(dp_s_star(:,cc),ctrs);
% %     h1 = h1./sum(h1);
% %     h2 = h2./sum(h2);
%     
%     subplot(nm,4,(cc-1)*4 + 4);
%     hold on
%     
%     plot(ctrs,h1,'Color',[0 0 0],'LineWidth',lw);
%     plot(ctrs,h2,'--','Color',[0.5 0.5 0.5],'LineWidth',lw);
%     
%     xlabel('stim. separability (min d'')','FontSize',fontSize_label);
%     ylabel('n. sims','FontSize',fontSize_label);
%     xlim([ctrs(1) ctrs(end)]);
%     set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
%     hold off
    
    %ctr_min = min([dp_s(:);dp_s_star(:)]);
    ctr_min = prctile([dp_s(:);dp_s_star(:)],0.25);
    ctr_max = max([dp_s(:);dp_s_star(:)]);
    ctrs = logspace(log10(ctr_min),log10(ctr_max), 31);
    h1 = hist(dp_s(:,cc),ctrs);
    h2 = hist(dp_s_star(:,cc),ctrs);
%     h1 = h1./sum(h1);
%     h2 = h2./sum(h2);
    
    subplot(nm,5,(cc-1)*5 + 4);
    
    semilogx(ctrs,h1,'Color',[0 0 0],'LineWidth',lw);
    hold on
    semilogx(ctrs,h2,'--','Color',[0.5 0.5 0.5],'LineWidth',lw);
    
    xlabel('stim. separability (min d'')','FontSize',fontSize_label);
    ylabel('n. sims','FontSize',fontSize_label);
    xlim([ctrs(1) ctrs(end)]);
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    hold off
    
    
    
    %%
    
    subplot(nm,5,(cc-1)*5 + 5);
    hold on

%     boxplot(perVar(:,:,cc)*100     ,{'time','stimulus','interaction'},'Positions',[1 4 7],'Width',0.75,'colors',[0 0 0]      ,'symbol','');
%     boxplot(perVar_star(:,:,cc)*100,{'time','stimulus','interaction'},'Positions',[2 5 8],'Width',0.75,'colors',[0.6 0.6 0.6],'symbol','');

    perI = 50;
    ps      = prctile(perVar(:,:,cc)*100     ,[50-perI/2 50+perI/2],1);
    ps_star = prctile(perVar_star(:,:,cc)*100,[50-perI/2 50+perI/2],1);
    
    ps_mu      = mean(perVar(:,:,cc)*100,1);
    ps_mu_star = mean(perVar_star(:,:,cc)*100,1);
    
    plot([1 4 7],ps_mu,'o','MarkerFaceColor',[0 0 0],'Color',[0 0 0])
    plot(ones(2,1)*[1 4 7],ps   ,'-','MarkerFaceColor',[0 0 0],'Color',[0 0 0])
    
    plot([1 4 7]+0.5,ps_mu_star,'o','MarkerFaceColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7]);
    plot(ones(2,1)*[1 4 7]+0.5,ps_star   ,'-','MarkerFaceColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7]);

    
    mm = ceil(max([perVar(:);perVar_star(:)])*20)/20 * 100;
    ylim([0 mm]);
    
    set(gca,'XTick',[1.5 4.5 7.5],'XTickLabel',{'time','stimulus','interaction'});
    
    xlabel('first component','FontSize',fontSize_label);
    ylabel('explained variance (%)','FontSize',fontSize_label);
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);
    hold off
end


set(gcf,'PaperSize',[15 (nm)*3],'PaperPosition',[0 0 15 (nm)*3]);
pause(1);
drawnow;
pause(1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResultsSummary.ss.mu      = mean(dp_s,1);
ResultsSummary.ss.mu_star = mean(dp_s_star,1);
ResultsSummary.ss.sig      = std(dp_s,1);
ResultsSummary.ss.sig_star = std(dp_s_star,1);


ResultsSummary.tt.mu      = mean(CC_t,1);
ResultsSummary.tt.mu_star = mean(CC_t_star,1);
ResultsSummary.tt.sig      = std(CC_t,1);
ResultsSummary.tt.sig_star = std(CC_t_star,1);


ResultsSummary.example.tt = nan(1,nm);
ResultsSummary.example.tt_star = nan(1,nm);

ResultsSummary.example.ss = nan(1,nm);
ResultsSummary.example.ss_star = nan(1,nm);
ResultsSummary.perVar = perVar;
ResultsSummary.perVar_star = perVar_star;


for ii = 1:nm
    [~,aa] = corrcoef(reshape(Ts(:,:,ii),[],1),tt_c(:));
    ResultsSummary.example.tt(ii) = aa(1,2)^2;
    [~,aa] = corrcoef([reshape(Ts(:,:,ii),[],1);reshape(Ts_star(:,:,ii),[],1)],[tt_c(:);tt_c_star(:)]);
    ResultsSummary.example.tt_star(ii) = aa(1,2)^2;
    
    dps = nan(C,C+C_star);
        
    for jj = 1:C
        A = rs{ii}.X_d((1:T)+(jj-1)*T,1,2);
        for kk = (jj+1):C
            B = rs{ii}.X_d((1:T)+(kk-1)*T,1,2);
            dps(jj,kk) = dPrime(A,B);
        end
        for kk = 1:C_star
            B = rs{ii}.X_star_d((1:T)+(kk-1)*T,1,2);
            dps(jj,C+kk) = dPrime(A,B);
        end
    end
    ResultsSummary.example.ss(ii)      = nanmin(nanmin(abs(dps(1:C,1:C))));
    ResultsSummary.example.ss_star(ii) = nanmin(nanmin(abs(dps(:,:))));
end

ResultsSummary.N = N_sims;


end


function [dp] = dPrime(A,B)
    A = A(:);
    B = B(:);
    dp = (mean(A)-mean(B))./sqrt(1/2*(var(A) + var(B)));
end



