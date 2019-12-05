

figDir = 'tex/figs_src';
mkdir(figDir);

saveResults = false;
useFigureComposer =  true;
if(useFigureComposer)
    addpath ~/FigureComposer/matlab %saves our figures using FigureComposer/DataNav format
end

% stimPoints = [1 5;
%               8 -2;
%               -3 -6];
stimPoints = [0  4;
              0   0;
              0 -4];
          
stimPoints(:,1) = 1;
          

theta =-20;
Mr = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
stimPoints = (Mr*stimPoints')';

stimPoints = stimPoints-mean(stimPoints,1);
          

timePoints = [-4 0;
              -2 0;
              0 0;
              2  0;
              4 0];
timePoints(:,2) = 0;

          
timePoints = timePoints-mean(timePoints,1);
          
rng(12122018);
NN = randn(15,2,2)*0.25;

% NN(1:5,:,2) = [-5 0;
%                -3 0;
%                 0 0;
%                 3 0;
%                 5 0];
% NN(6:10,:,2) = [1 -6;
%                0.5 -3;
%                 0 0;
%                 -0.5 3;
%                 -1 6];0
% NN(11:15,:,2) = [-0.5 2;
%                -0.25 1;
%                 0 0;
%                 0.25 -1;
%                 0.5 -2];
NN(1:5  ,:,2) = [-2    -3;
                 -1   -1.5;
                 0     0;
                  1   1.5;
                  2   3];
NN(11:15 ,:,2) = [ 2    2;
                   1  1;
                  0    0 ;
                  -1  -1;
                  -2  -2];
NN(6:10 ,:,2) = -(NN(1:5  ,:,2) + NN(11:15,:,2));
            
            
% NN(1:5,:,2) = -[5 -5;
%                3 -3;
%                0 0;
%                -3 3;
%                -5 5];
% NN(6:10,1,2) = -NN(1:5,1,2);
% NN(6:10,2,2) = -NN(1:5,2,2);
% NN(11:15,:,2) = [0 0;
%                0 0;
%                 0 0;
%                 0.0 0;
%                 0.0 0];

ccs = [1 0 0;
       0 1 0;
       0 0 1];
          
ms = 5;
lw  = 1.5;
lw2 = 1;

fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;


NR = 4;
NC = 4;
          
figure(1);
clf
for ff = 1:2
    subplot(NR,NC,1+(ff-1)*NC*2)
    hold on
    for ii = 1:3
        plot(stimPoints(ii,1),stimPoints(ii,2),'o','MarkerFaceColor',ccs(ii,:),'Color',ccs(ii,:),'MarkerSize',ms);
    end

    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('stimulus component','fontSize',fontSize_title);
    axis square
    hold off


    subplot(NR,NC,2+(ff-1)*NC*2)
    hold on
    plot(timePoints(:,1),timePoints(:,2),'o-','MarkerFaceColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7],'MarkerSize',ms);
    plot(timePoints(1,1),timePoints(1,2),'o-','Color',[0.7 0.7 0.7]*0,'MarkerSize',ms);

    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('time component','fontSize',fontSize_title);
    axis square
    hold off

    subplot(NR,NC,3+(ff-1)*NC*2);
    hold on

    for ii = 1:3
        tts  = (1:5) + (ii-1)*5;
        plot(NN(tts,1,ff),NN(tts,2,ff),'o-','MarkerFaceColor',ccs(ii,:),'Color',ccs(ii,:),'MarkerSize',ms);
        plot(NN(tts(1),1,ff),NN(tts(1),2,ff),'o-','MarkerFaceColor',ccs(ii,:),'Color',[0 0 0],'MarkerSize',ms);
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('interaction','fontSize',fontSize_title);
    axis square
    hold off

    subplot(NR,NC,4+(ff-1)*NC*2)
    hold on

    pps = nan(5,2,3);
    
    for ii = 1:3
        tts  = (1:5) + (ii-1)*5;
        timePoints_c = timePoints;
        
        pps(:,:,ii) = [NN(tts,1,ff)+stimPoints(ii,1) + timePoints_c(:,1), NN(tts,2,ff)+stimPoints(ii,2) + timePoints_c(:,2)]; 
        plot(pps(:,1,ii),pps(:,2,ii),'o-','MarkerFaceColor',ccs(ii,:),'Color',ccs(ii,:),'MarkerSize',ms);
        plot(pps(1,1,ii),pps(1,2,ii),'o-','MarkerFaceColor',ccs(ii,:),'Color',[0 0 0],'MarkerSize',ms);
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('total','fontSize',fontSize_title);
    axis square
    hold off
    
    
    subplot(NR,NC,5+(ff-1)*NC*2)
    hold on

    
    for ii = 1:3
        plot(mean(pps(:,1,ii)),mean(pps(:,2,ii)),'o','MarkerFaceColor',ccs(ii,:),'Color',ccs(ii,:),'MarkerSize',ms);
    end
    
    
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('stimulus mean','fontSize',fontSize_title);
    axis square
    hold off

    subplot(NR,NC,6+(ff-1)*NC*2)
    hold on

    
    plot(mean(pps(:,1,:),3),mean(pps(:,2,:),3),'o-','MarkerFaceColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7],'MarkerSize',ms);

    
    
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('time mean','fontSize',fontSize_title);
    axis square
    hold off
end
plotSize = 2;
set(gcf,'PaperSize',[NC NR]*plotSize,'PaperPosition',[0 0 NC NR]*plotSize);
if(saveResults)
    if(useFigureComposer)
        matfig2fyp(gcf,sprintf('%s/cartoon_rough.fyp',figDir));
    else
        saveas(gcf,sprintf('%s/cartoon_rough.eps',figDir),'epsc');
    end
end
%%

X_gamma = nan(15,2,3);
X_gamma(:,:,1) = [repmat(stimPoints(1,:),5,1);
                  repmat(stimPoints(2,:),5,1);
                  repmat(stimPoints(3,:),5,1)];
X_gamma(:,:,2) = repmat(timePoints,3,1);
X_gamma(:,:,3) = NN(:,:,2);

Xs = sum(X_gamma,3);

lambda_k = 1e-6;
lengthScale = 5;

aa = -10:0.1:10;
[xx,yy] = meshgrid(aa,aa);
Xs_star = [xx(:),yy(:)];

[results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA_sqExp(Xs,X_gamma,lambda_k,lengthScale,Xs_star,1);
%[results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d] = kdPCA_md(Xs,X_gamma,lambda_k,lengthScale,Xs_star,1);

figure(2);
clf
for ii = 1:2
    subplot(1,2,ii);
    C = reshape(results_kdPCA.X_star_d(:,1,ii),[],length(aa));
    
    C = max(C,-5);
    C = min(C, 5);
    hold on
    imagesc(aa,aa,C);
    colorbar;
    
    
    for jj = 1:3
        tts  = (1:5) + (jj-1)*5;
        
        pps = Xs(tts,:); 
        plot(pps(:,1),pps(:,2),'o-','MarkerFaceColor',ccs(jj,:),'Color',ccs(jj,:),'MarkerSize',ms);
        plot(pps(1,1),pps(1,2),'o-','MarkerFaceColor',ccs(jj,:),'Color',[0 0 0],'MarkerSize',ms);
    end
    hold off
    xlim([-10 10]);
    ylim([-10 10]);
    
end