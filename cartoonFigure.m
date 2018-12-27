

figDir = 'tex/figs_src';
mkdir(figDir);

saveResults = false;
useFigureComposer =  false;
if(useFigureComposer)
    addpath ~/FigureComposer/matlab %saves our figures using FigureComposer/DataNav format
end

stimPoints = [1 5;
              8 2;
              -4 -7];

timePoints = [-5 1;
              -1 0.1;
              -0.5 -1;
              0   -0.1;
              0.75 0.5];
rng(12122018);
NN = randn(15,2)*0.25;

ccs = [1 0 0;
       0 1 0;
       0 0 1];
          
ms = 5;
lw  = 1.5;
lw2 = 1;

fontSize_axis = 10;
fontSize_label = 12;
fontSize_title = 12;


MM = zeros(2,2,3);
MM(:,:,1) = [1 0;
             0 1];
MM(:,:,2) = [2.5 0;
             0    0.5];
MM(:,:,3) = [0.5 0;
             0    0.5];
R1 = deg2rad(45);
R2 = deg2rad( -45);
         
         
MM(:,:,2) = [cos(R1) -sin(R1);  
             sin(R1)  cos(R1)]*MM(:,:,2);
MM(:,:,3) = [cos(R2) -sin(R2);  
             sin(R2)  cos(R2)]*MM(:,:,3);
          
figure(1);
clf
for ff = 1:2
    subplot(2,4,1+(ff-1)*4)
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


    subplot(2,4,2+(ff-1)*4)
    hold on

    if(ff == 2)
        for ii = 1:3
            timePoints_c = timePoints;
            timePoints_c = (MM(:,:,ii)*timePoints_c')';

            plot(timePoints_c(:,1),timePoints_c(:,2),'o-','MarkerFaceColor',ccs(ii,:),'Color',[0.7 0.7 0.7],'MarkerSize',ms);
        end
    else
        plot(timePoints(:,1),timePoints(:,2),'o-','MarkerFaceColor',[0.7 0.7 0.7],'Color',[0.7 0.7 0.7],'MarkerSize',ms);
    end

    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('time component','fontSize',fontSize_title);
    axis square
    hold off

    subplot(2,4,3+(ff-1)*4);
    hold on

    for ii = 1:3
        tts  = (1:5) + (ii-1)*5;
        plot(NN(tts,1),NN(tts,2),'o','MarkerFaceColor',ccs(ii,:),'Color',ccs(ii,:),'MarkerSize',ms);
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('noise','fontSize',fontSize_title);
    axis square
    hold off

    subplot(2,4,4+(ff-1)*4)
    hold on

    for ii = 1:3
        tts  = (1:5) + (ii-1)*5;
        timePoints_c = timePoints;
        if(ff == 2)
            timePoints_c = (MM(:,:,ii)*timePoints_c')';
        end
        
        plot(NN(tts,1)+stimPoints(ii,1) + timePoints_c(:,1),NN(tts,2)+stimPoints(ii,2) + timePoints_c(:,2),'o-','MarkerFaceColor',ccs(ii,:),'Color',ccs(ii,:),'MarkerSize',ms);
    end
    set(gca,'TickDir','out','box','off','FontSize',fontSize_axis);

    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('latent dim. 1','FontSize',fontSize_label);
    ylabel('latent dim. 2','FontSize',fontSize_label);
    title('total','fontSize',fontSize_title);
    axis square
    hold off

    set(gcf,'PaperSize',[8 4],'PaperPosition',[0 0 8 4]);
end
if(saveResults)
    if(useFigureComposer)
        matfig2fyp(gcf,sprintf('%d/cartoon_rough.fyp',figDir));
    else
        saveas(gcf,sprintf('%d/cartoon_rough.eps',figDir),'epsc');
    end
end
          