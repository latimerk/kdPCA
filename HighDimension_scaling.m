rng(20190512);
saveFig = false;

D =  6;
N = 50;
T_0 = 10;
T = T_0*D;
NE = 8;

ms = 6;
lw  = 1.5;
lw2 = 1;

lambda = 1;
lambda_k = 1;
lengthScale = 5;
R = 4;

fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;

sig = 1;

X_0 = zeros(D*T_0,D);

for dd = 1:D
    for tt = 1:T_0
        ii = tt + (dd-1)*T_0;
        
        X_0(ii,dd) = tt;
        X_0(ii,1:(dd-1)) = T_0;
        
        
%         X_0(ii,1:(dd-2)) = 0;
%         if(dd > 1)
%             X_0(ii,(dd-1)) = max(0,T_0-tt);
%         end
%         
%         
%         X_0(ii,1:(dd-3)) = 0;
%         if(dd > 1)
%             X_0(ii,(dd-1)) = T_0;
%         end
%         if(dd > 2)
%             X_0(ii,(dd-2)) = max(0,T_0-tt);
%         end
    end
end
X_0 = X_0 - ones(1,D)*T_0/2;


      
scales = [linspace(0.5 ,2   ,D);  
          linspace(0.75,1.5 ,D);  
          linspace(1   ,1   ,D);  
          linspace(1.5 ,0.75,D);  
          linspace(2   ,0.5 ,D)];
      
      
% scales = repmat([0.5 0.75 1 1.5 2]',[1 D]);
      
Ls = zeros(3*D*T_0,D);

Ls((1:(T_0*D)) + (T_0*D)*0,:) = X_0.*scales(1,:);
Ls((1:(T_0*D)) + (T_0*D)*1,:) = X_0.*scales(3,:);
Ls((1:(T_0*D)) + (T_0*D)*2,:) = X_0.*scales(5,:);
      
Ls_star = zeros(2*D*T_0,D);

Ls_star((1:(T_0*D)) + (T_0*D)*0,:) = X_0.*scales(2,:);
Ls_star((1:(T_0*D)) + (T_0*D)*1,:) = X_0.*scales(4,:);
      
C_star = 2;
C      = 3;

if(D == N)
    W = eye(N);
else
    W = randn(D,N);
end
Xs_0 = randn(C*T,N)*sig+Ls*W;
Xs_star_0 = randn(C_star*T,N)*sig + Ls_star*W;

[Ys,mu,st] = zscore(Xs_0);
st(st==0) = 1;
Ys_star = (Xs_star_0-mu)./st;

Xs = Ys;
Xs_star = Ys_star;
N2 = size(Xs,2);

%%


X_gamma = zeros(T*C,N2,3);

for ii = 1:C
    tts_ii = (1:T)+(ii-1)*T;
    X_gamma(tts_ii,:,2) = ones(length(tts_ii),1)*mean(Xs(tts_ii,:),1);
end
for ii = 1:T
    b1 = ii:T:(C*T);
    X_gamma(b1,:,1) = ones(length(b1),1)*mean(Xs(b1,:),1);
end

X_gamma(:,:,3) = Xs-sum(X_gamma(:,:,1:2),3);





[COEFF1, SCORE1, LATENT1, TSQUARED1, EXPLAINED1,MU1] = pca(Xs);

SCORE1_star = (Xs_star-MU1)*COEFF1;


[results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA_sqExp(Xs,X_gamma,lambda_k,lengthScale,Xs_star,R);
[ results_dPCA.X_d, results_dPCA.X_star_d, results_dPCA.X_hat_d, results_dPCA.X_star_hat_d] = dPCA(Xs,X_gamma,lambda,Xs_star,R);

figure(1);
clf
sets = [1 2;
        1 3;
        2 3];
for gg = 1:size(sets,1)
    subplot(3,4,gg);
    hold on
    L = SCORE1(:,sets(gg,:));
    L_star = SCORE1_star(:,sets(gg,:));
    min_c = min(min(min(L_star(:,1:2))),min(min(L(:,1:2))));
    max_c = max(max(max(L_star(:,1:2))),max(max(L(:,1:2))));

    colors_c = zeros(C,3);
    colors_c_star = zeros(C_star,3);
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        pp = plot(L(tts,1),L(tts,2),'-','LineWidth',lw);
        colors_c(ii,:) = pp.Color;
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        pp = plot(L_star(tts,1),L_star(tts,2),'--','LineWidth',lw);
        colors_c_star(ii,:) = pp.Color;
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        pp = plot(L_star(tts(1),1),L_star(tts(1),2),'o','LineWidth',lw2,'Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'MarkerSize',ms);
    end
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        pp = plot(L(tts(1),1),L(tts(1),2),'o','LineWidth',lw2,'Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'MarkerSize',ms);
    end
    
    lims = [floor(min_c*2)./2 ceil(max_c*2)./2];
    xlim(lims);
    ylim(lims);
    set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
    
    xlabel(sprintf('PC %d',sets(gg,1)),'FontSize',fontSize_label);
    ylabel(sprintf('PC %d',sets(gg,2)),'FontSize',fontSize_label);
    
    hold off
end
subplot(3,4,4)
plot(1:min(NE,length(EXPLAINED1)),EXPLAINED1(1:min(NE,length(EXPLAINED1))),'o-');
xlabel('PC','FontSize',fontSize_label);
ylabel('% Variance Explained','FontSize',fontSize_label);
set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
xlim([1 min(NE,length(EXPLAINED1))]);

lims = [-7.5  8.5;
        -4.5 4;
        -3   2];

for gg = 1:3
    %%
    dds = [1 2];
    subplot(3,4,4+gg);
    hold on
    L = results_dPCA.X_d(:,dds,gg);
    L_star = results_dPCA.X_star_d(:,dds,gg);
    
    min_c = min(min(min(L_star(:,1:2))),min(min(L(:,1:2))));
    max_c = max(max(max(L_star(:,1:2))),max(max(L(:,1:2))));
    fprintf('gg = %d\n',gg)
    fprintf('   min = %.2f\n',min_c);
    fprintf('   max = %.2f\n',max_c);

    colors_c = zeros(C,3);
    colors_c_star = zeros(C_star,3);
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        pp = plot(L(tts,1),L(tts,2),'-','LineWidth',lw);
        colors_c(ii,:) = pp.Color;
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        pp = plot(L_star(tts,1),L_star(tts,2),'--','LineWidth',lw);
        colors_c_star(ii,:) = pp.Color;
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        pp = plot(L_star(tts(1),1),L_star(tts(1),2),'o','LineWidth',lw2,'Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'MarkerSize',ms);
    end
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        pp = plot(L(tts(1),1),L(tts(1),2),'o','LineWidth',lw2,'Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'MarkerSize',ms);
    end
    set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
    
    lims(gg,:) = [floor(min_c*2)./2 ceil(max_c*2)./2];
    
    xlim(lims(gg,:));
    ylim(lims(gg,:));
    
    hold off
end

lims = [-7 8;
        -3 3;
        -4.5 4];

for gg = 1:3
    
    %%
    subplot(3,4,8+gg);
    hold on
    
    dds = [1 2];
    
    L = results_kdPCA.X_d(:,dds,gg);
    L_star = results_kdPCA.X_star_d(:,dds,gg);
    colors_c = zeros(C,3);
    colors_c_star = zeros(C_star,3);
    
    
    min_c = min(min(min(L_star(:,1:2))),min(min(L(:,1:2))));
    max_c = max(max(max(L_star(:,1:2))),max(max(L(:,1:2))));
    fprintf('gg = %d\n',gg)
    fprintf('   min = %.2f\n',min_c);
    fprintf('   max = %.2f\n',max_c);
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        pp = plot(L(tts,1),L(tts,2),'-','LineWidth',lw);
        colors_c(ii,:) = pp.Color;
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        pp = plot(L_star(tts,1),L_star(tts,2),'--','LineWidth',lw);
        colors_c_star(ii,:) = pp.Color;
    end
    for ii = 1:C_star
        tts = (1:T) + (ii-1)*T;
        pp = plot(L_star(tts(1),1),L_star(tts(1),2),'o','LineWidth',lw2,'Color',colors_c_star(ii,:),'MarkerFaceColor',[0 0 0],'MarkerSize',ms);
    end
    for ii = 1:C
        tts = (1:T) + (ii-1)*T;
        pp = plot(L(tts(1),1),L(tts(1),2),'o','LineWidth',lw2,'Color',colors_c(ii,:),'MarkerFaceColor',[0 0 0],'MarkerSize',ms);
    end
    set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
    lims(gg,:) = [floor(min_c*2)./2 ceil(max_c*2)./2];
    xlim(lims(gg,:));
    ylim(lims(gg,:));
    hold off
end

if(saveFig)
    set(gcf,'PaperSize',[12 9],'PaperPosition',[0 0 12 9]);
    addpath ~/FigureComposer/
    addpath ~/FigureComposer/matlab/
    figDir = 'tex/figs_src';
    matfig2fyp(gcf,sprintf('%s/HighDsim_raw.fyp',figDir));
end
