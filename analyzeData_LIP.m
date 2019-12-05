
load('Data/lipProcessed_marginals.mat'); %if file does not exist, run getLIPmarginals.m


useFigureComposer = true;
saveResults = false;
figDir = 'tex/figs_src';

if(useFigureComposer)
    addpath ~/FigureComposer/matlab/
end

lengthScale = 50;
R = 5;

ms = 6;

fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;


lambda_k = 1;

[results_kdPCA.X_d,results_kdPCA.X_star_d,results_kdPCA.X_hat_d,results_kdPCA.X_star_hat_d,~,~,results_kdPCA.hhs] = kdPCA_sqExp(X,X_gamma,lambda_k,lengthScale,X_star,R);


lambda   = 1;
[ results_dPCA.X_d, results_dPCA.X_star_d, results_dPCA.X_hat_d, results_dPCA.X_star_hat_d,results_dPCA.dds,results_dPCA.ffs] = dPCA(X,X_gamma,lambda,X_star,R);

%%
lambda2s = [0.1 10 100];

results_dPCA2 = cell(1,length(lambda2s));

for ii = 1:length(lambda2s)
    results_dPCAc = struct();
    [ results_dPCAc.X_d, results_dPCAc.X_star_d, results_dPCAc.X_hat_d, results_dPCAc.X_star_hat_d,results_dPCAc.dds,results_dPCAc.ffs] = dPCA(X,X_gamma,lambda2s(ii),X_star,R);
    results_dPCA2{ii} = results_dPCAc;
end

%%
lw1 = 1;
lw2 = 1;

tt_axis = (1:T)*0.001 - 0.2;



xlims = [-90 55;
        -30 30;
        -72 68;
        -12 12];

ylims = [-55 90;
        -30 30;
        -72 68;
        -12 12];

xticks = {-80:40:40;
          -30:15:30;
          -60:30:60;
          -12:6:12};
yticks = {-40:40:80;
          -30:15:30;
          -60:30:60;
          -12:6:12};

NC = 6;
%ls_test = {'-.'           ,':'             ,'--'          ,  '-'};
%          {Out-RF,T-Flash}, {In-RF,T-Flash}, {Out-RF,T-ON}, {In-RF, T-ON}
%          [1 1]              [2 1]           [1 2]           [2 2]
lss = {'-.','--';
       ':','-'};
colors = [ 255 0 0;
      180 0 0;
      0 0 0;
      0  0 180;
      0  0 255]./255;

colors_test = [colors(Test_Cs(1),:);
               colors(Test_Cs(2),:);
               colors(Test_Cs(1),:);
               colors(Test_Cs(2),:)];
strs2 = {'time','decision','targets','interaction'};
ls_test = {'-.',':','--','-'};
                                   
for pp = 1:2
    
    if(pp == 1)
        rs = {results_dPCA,results_kdPCA};
        strs = {'dPCA','kdPCA'};

    else
        rs = results_dPCA2;
        strs = cell(length(rs),1);
        for ii = 1:length(rs)
            strs{ii} = sprintf('dPCA \\lambda=%.1f',lambda2s(ii));
        end
    end
    RR = length(rs);

    figure(20+pp);
    clf

    lims = [-90 50;
            -20 25;
            -22 22;
            -10 50];
        
    for dd = 1:RR
        for cc = 1:4
            subplot(RR,NC,(dd-1)*NC+cc)
            hold on
            if(lw1 > 0)
                for ii = 1:C
                    for jj = 1:D
                        for kk = 1:E
                            plot(rs{dd}.X_d(idxs{ii,jj,kk},1,cc),rs{dd}.X_d(idxs{ii,jj,kk},2,cc),lss{jj,kk},'Color',colors(Train_Cs(ii),:),'LineWidth',lw1)
                        end
                    end
                end
            end
            
            if(lw2 > 0)
                for hh = 1:H
                    plot(rs{dd}.X_star_d((1:T)+(hh-1)*T,1,cc),rs{dd}.X_star_d((1:T)+(hh-1)*T,2,cc),ls_test{hh},'Color',colors_test(hh,:),'LineWidth',lw2);
                end
            end
            
            set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
            %title(sprintf('%s space: %s',strs2{cc},strs{dd}),'FontSize',fontSize_title);
            if(dd == RR)
                xlabel('component 1','FontSize',fontSize_label);
            end
            if(cc == 1)
                ylabel('component 2','FontSize',fontSize_label);
            end
            axis square;
            axis equal;
            
            
            xlim(xlims(cc,:));
            ylim(ylims(cc,:));
            set(gca,'XTick',xticks{cc},'YTick',yticks{cc});
            hold off


        end
        for ll = 1:2
            subplot(RR,NC,(dd-1)*NC + 4 + ll);
            hold on
            
            for ii = 1:C
                for jj = 1:D
                    for kk = 1:E
                        plot(tt_axis,rs{dd}.X_d(idxs{ii,jj,kk},ll,cc),lss{jj,kk},'Color',colors(Train_Cs(ii),:),'LineWidth',lw1);
                    end
                end
            end
            
            if(lw2 > 0)
                for hh = 1:H
                    plot(tt_axis,rs{dd}.X_star_d((1:T)+(hh-1)*T,ll,cc),ls_test{hh},'Color',colors_test(hh,:),'LineWidth',lw2);
                end
            end

            lls = [-12 12];
            %plot([0 0],lls,'-','LineWidth',0.5,'Color',[0.7 0.7 0.7]);
            set(gca,'TickDir','out','box','off','fontSize',fontSize_axis);
            if(dd == RR)
                xlabel('time from motion onset (s)','FontSize',fontSize_label);
            end
            ylabel(sprintf('component %d',ll),'FontSize',fontSize_label);
            xlim([-0.2 0.5]);
            ylim(lls);

            set(gca,'YTick',-12:6:12,'XTick',-0.2:0.2:0.4);
            hold off
        end
    end

    drawnow;

    set(gcf,'PaperSize',[6 RR]*3,'PaperPosition',[0 0 6 RR]*3);

    if(saveResults)
        if(useFigureComposer)
            matfig2fyp(gcf,sprintf('%s/DataLIP_%d_rough.fyp',figDir,pp));
        else
            saveas(gcf,sprintf('%s/DataLIP_1_rough.eps',figDir),'epsc');
        end
    end
end