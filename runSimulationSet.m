
N = 50;

saveResults = true;
useFigureComposer = true;
figDir = 'tex/figs_src';
mkdir(figDir);

if(useFigureComposer)
    addpath ~/FigureComposer/matlab %saves our figures using FigureComposer/DataNav format
end

rng(12122018);
figs_plainSimulation;
if(saveResults)
    if(useFigureComposer)
        drawnow;
        figure(1);
        matfig2fyp(gcf,sprintf('%s/plainSim_1_raw.fyp',figDir));
        figure(2);
        matfig2fyp(gcf,sprintf('%s/plainSim_2_raw.fyp',figDir));
    else
        figure(1);
        saveas(gcf,sprintf('%s/plainSim_1_raw.eps',figDir),'epsc');
        figure(2);
        saveas(gcf,sprintf('%s/plainSim_2_raw.eps',figDir),'epsc');
    end
    save('ResultsSummary_plain.mat','ResultsSummary_plain');
end

%%
rng(12132018);

figs_scaledSimulation;

if(saveResults)
    if(useFigureComposer)
        drawnow;
        figure(1);
        matfig2fyp(gcf,sprintf('%s/scaledSim_1_raw.fyp',figDir));
        figure(2);
        matfig2fyp(gcf,sprintf('%s/scaledSim_2_raw.fyp',figDir));
        figure(3);
        matfig2fyp(gcf,sprintf('%s/scaledSim_3_raw.fyp',figDir));
    else
        drawnow;
        figure(1);
        saveas(gcf,sprintf('%s/scaledSim_1_raw.eps',figDir),'epsc');
        figure(2);
        saveas(gcf,sprintf('%s/scaledSim_2_raw.eps',figDir),'epsc');
        figure(3);
        saveas(gcf,sprintf('%s/scaledSim_3_raw.eps',figDir),'epsc');
    end

    save('ResultsSummary_scaled.mat','ResultsSummary_scaled');
end
%%
rng(04282019); 

figs_rotatedSimulation;

if(saveResults)
    if(useFigureComposer)
        drawnow;
        figure(1);
        matfig2fyp(gcf,sprintf('%s/rotatedSim_1_raw.fyp',figDir));
        figure(2);
        matfig2fyp(gcf,sprintf('%s/rotatedSim_2_raw.fyp',figDir));
        figure(3);
        matfig2fyp(gcf,sprintf('%s/rotatedSim_3_raw.fyp',figDir));
    else
        drawnow;
        figure(1);
        saveas(gcf,sprintf('%s/rotatedSim_1_raw.eps',figDir),'epsc');
        figure(2);
        saveas(gcf,sprintf('%s/rotatedSim_2_raw.eps',figDir),'epsc');
        figure(3);
        saveas(gcf,sprintf('%s/rotatedSim_3_raw.eps',figDir),'epsc');
    end

    save('ResultsSummary_rotated.mat','ResultsSummary_rotated');
end