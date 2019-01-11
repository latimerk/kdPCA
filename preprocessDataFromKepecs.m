% Data processing for kdPCA example in
%   Nonlinear demixed component analysis for neural population data as a low-rank kernel regression problem
%
%   by Kenneth W. Latimer (2019)
%
% This is a slightly modified script from
%
%   preprocess_ClaudiaAdam.m
%
%   D Kobak+, W Brendel, C Constantinidis, CE Feierstein, A Kepecs, ZF Mainen, X-L Qi, R Romo, N Uchida, CK Machens
%   Demixed principal component analysis of neural population data
%   eLife, 2016
%
%   Downloaded in Dec 2018 from https://github.com/machenslab/elife2016dpca 
%   
%

clear all

% datadir = '/home/dmitry/Dropbox/Machens Lab Sharing/Data/ClaudiaFeierstein/cellbase1/';
%datadir = '/home/dmitry/Dropbox/Machens Lab Sharing/Data/Adam Kepecs/';
datadir = '~/gitCode/kdPCA/Data/data-Kepecs/'; %directory containing the ofc-2 data (subfolders are N1, N48, and N49)

% outputFileName = 'data_claudia_piecewise.mat';
outputFileName = 'data_adam.mat';



restretch = true;
piecewiseStitch = ~restretch;

%%

alignment = {'OdorPokeIn','OdorPokeOut','WaterPokeIn','WaterValveOn','WaterPokeOut'};
tbefore = 0.5;       % start this many secs before first alignment variable
tafter = 0.5;        % end this many secs after the last alignment variable

% get all sessions' paths (each session is a folder)
sessions = {};
sessionsRats = [];
dirs = dir(datadir);    % note that first two are '.' and '..'
for i=3:length(dirs)    % loop over rats
    subdirs = dir([datadir dirs(i).name]);
    for j=3:length(subdirs)
        if subdirs(j).isdir
            sessions{end+1} = [datadir dirs(i).name '/' subdirs(j).name '/'];
            sessionsRats = [sessionsRats i-2];
        end
    end
end
display(['Found ' num2str(length(sessions)) ' sessions from ' num2str(length(unique(sessionsRats))) ' rats'])

% compute median alignment timepoints
numOfTrialsUntilNow = 0;
alignmenttimes = [];
odoursTable = [];
neuronsPerSession = zeros(length(sessions),1);
for s = 1:length(sessions)
    load([sessions{s}, 'TrialEvents2.mat']);
    
    if length(WaterPokeIn) ~= length(WaterValveOff) % skips two broken sessions
        display(['Skipping session ' sessions{s}])
        continue
    end
    
    if isnan(sum(OdorCategory(~isnan(OdorValveID)))) % skips a session with 50/50 mixtures
        display(['Skipping session ' sessions{s}])
        continue
    end
    
    numOfTrials = length(OdorPokeValid);
    for a = 1:length(alignment)
        alignmenttimes(numOfTrialsUntilNow+1:numOfTrialsUntilNow+numOfTrials, a) = eval(alignment{a})';
    end
    
    if ~isempty(strfind(datadir, 'Feierstein'))
        ThisSessionDecisions = WaterPokeID;
    end
    if ~isempty(strfind(datadir, 'Kepecs'))
        ThisSessionDecisions = ChoiceDir;
    end
    
    % WaterValveOn -- 4
    ratCorrectAndGotReward = find(Correct == 1 & ~isnan(WaterValveOn) & ~isnan(ThisSessionDecisions));
    meanWaterValveOn(s) = mean(WaterValveOn(ratCorrectAndGotReward) - WaterPokeIn(ratCorrectAndGotReward));
    %alignmenttimes(numOfTrialsUntilNow + find(Error == 1), 4) = alignmenttimes(numOfTrialsUntilNow + find(Error == 1), 3) ...
    %                                                            + meanWaterValveOn(s);
    
    numOfTrialsUntilNow = numOfTrialsUntilNow + numOfTrials;
    
    % which odours were used?
    odours = unique(OdorValveID(~isnan(OdorValveID)));
    for o=1:length(odours)
        odoursTable(s, odours(o)) = OdorCategory(find(OdorValveID==odours(o),1));
    end
    
    % count neurons
    files = dir(sessions{s});
    filenames = {files.name};
    for f = 1:length(filenames)
        if ~isempty(regexp(filenames{f}, '^Sc.*\.mat$', 'once'))            
            neuronsPerSession(s) = neuronsPerSession(s) + 1;
        end
    end
    
    % if mixtures were used (only works for Claudia's data, not for Adam's)
    if exist('protocol', 'var')
        pr{s} = protocol;
    end
    
    if exist('OdorRatio', 'var')
        ratiosUsed(s, unique(OdorRatio)+1) = 1;
    end
end

display(['Found ' num2str(sum(neuronsPerSession)) ' neurons in total'])

figure
imagesc(odoursTable)
hold on
for r = find(diff(sessionsRats))
    plot([0.5 size(odoursTable, 2)+0.5], [r r]+0.5, 'w', 'LineWidth', 3)
end
xlabel('Odour id')
ylabel('Sessions')

if exist('ratiosUsed', 'var')
    figure
    imagesc(ratiosUsed)
    hold on
    for r = find(diff(sessionsRats))
        plot([0.5 100.5], [r r]+0.5, 'w', 'LineWidth', 3)
    end
    xlabel('Odour ratio used')
    ylabel('Sessions')
end

if exist('protocol', 'var')
    sessionsPure = ones(1,length(sessions));   % Clau's data
    sessionsPure(strcmp(pr,'binmix')) = 0;
else
    sessionsPure = zeros(1, length(sessions)); % Adam's data => all mixtures
end

% figure
% plot(neuronsPerSession)
% hold on
% xlabel('Session')
% ylabel('Number of neurons')

% selecting the trials with all events present
% display([num2str(length(find(~isnan(sum(alignmenttimes,2))))) ' fully finished trials out of ' num2str(size(alignmenttimes,1))])
alignmentMedian = nanmedian(alignmenttimes);
alignmentMedian = alignmentMedian - alignmentMedian(1);

% figure
% for a = 1:length(alignment)
%     subplot(length(alignment),1,a)
%     hold on
%     hist(alignmenttimes(:,a), 1000);
%     plot(median(alignmenttimes(:,a)), 0, 'r.', 'MarkerSize', 40)
%     plot(mean(alignmenttimes(:,a)), 0, 'g.', 'MarkerSize', 40)
%     xlim([-10 50]);
%     title(alignment{a})
% end

time = -tbefore:0.01:(alignmentMedian(end)+tafter);
% time = (-100:455)*0.01;
x = -50:50;
gaussKernel = 1/sqrt(2*pi)/5 * exp(-x.^2/5^2/2) * 10^2;    %5*10ms = 50ms width

unit = 0; % unit counter
unstableUnits = 0;

% loop over sessions
for s = 1:length(sessions)
    display(['Processing session ' num2str(s) ' out of ' num2str(length(sessions)) '...'])
    pause(0.001)
    
    % reload data for experimental parameters
    load([sessions{s} 'TrialEvents2.mat']);
    
    if length(WaterPokeIn) ~= length(WaterValveOff) % skips two broken sessions
        display('Skipping this session due to broken number of trials')
        continue
    end
    
    if isnan(sum(OdorCategory(~isnan(OdorValveID)))) % skips one broken session
        display('Skipping this session due to missing odour categories')
        continue
    end
    
    odours = unique(OdorValveID(~isnan(OdorValveID)));
    numOfTrials = length(OdorPokeValid);

    alignmentSession = [];
    for a = 1:length(alignment)
        alignmentSession(:, a) = eval(alignment{a})';
    end
    
    if ~isempty(strfind(datadir, 'Feierstein'))
        minAnticipationLength = 0.2;
    end
    if ~isempty(strfind(datadir, 'Kepecs'))
        minAnticipationLength = 0.3;
    end
    
    alignmentSession(find(Error == 1), 4) = alignmentSession(find(Error == 1), 3) + minAnticipationLength;
    
    % find correct trials in this session
    correctTrials = [];
    stimulus = [];
    stimulusCategory = [];
    stimulusRatio = [];
    decision = [];
    reward = [];
    for t = 1:numOfTrials
        % in Claudia's data rat #3 had a new pair of odours introduced
        % later in each session. Following Claudia's Neuron paper, we 
        % are disregarding these later parts of each session
        if ~isempty(strfind(datadir, 'Feierstein')) && sessionsRats(s)==3 ...
                && t == find(OdorValveID==2 | OdorValveID==6, 1)
            break
        end       
        
        if ~isempty(strfind(datadir, 'Feierstein'))
           ThisSessionDecisions = WaterPokeID;
        end
        if ~isempty(strfind(datadir, 'Kepecs'))
           ThisSessionDecisions = ChoiceDir;
        end
        
        [~, id] = sort(alignmentSession(t,:));
        if OdorPokeValid(t) && ~isnan(sum(alignmentSession(t,:))) ...
                && max(abs(id-(1:size(alignmentSession,2))))==0 && ...
                ~isnan(ThisSessionDecisions(t)) && ThisSessionDecisions(t)>0 && ~isnan(OdorValveID(t)) ...
                && ~(Correct(t)==1 && WaterValveID(t) == 0) % trial is correct, but WaterValveID==0, i.e. rat did not wait for the reward
            correctTrials(end+1) = t;
            stimulus(end+1) = OdorValveID(t);
            decision(end+1) = ThisSessionDecisions(t);   
            stimulusCategory(end+1) = OdorCategory(t);
            if ThisSessionDecisions(t) == OdorCategory(t)
                assert(Correct(t)==1)
                reward(end+1) = 1;
            else
                assert(Error(t)==1)
                reward(end+1) = 2;      % 2 means no reward (2~0 mod. 2)
            end
            
            if ~isempty(strfind(datadir, 'Kepecs'))
                ratiosAvailable = [0    32    44    56    68    100];
                stimulusRatio(end+1) = find(ratiosAvailable == OdorRatio(t));
            end
        end
    end
            
    % loop over units
    files = dir(sessions{s});
    filenames = {files.name};
    
    for f = 1:length(filenames)
        if isempty(regexp(filenames{f}, '^Sc.*\.mat$', 'once'))
            % not a spiketrain
            continue
        end
        
        load([sessions{s} filenames{f}])
        spiketrain = TS/10000; % convert to seconds
        
        % Checking if this neuron is firing at all
        trialSpikes = zeros(1, length(correctTrials));
        for tr = 1:length(correctTrials)
            sb = TrialStart(correctTrials(tr)) + alignmentSession(correctTrials(tr), 1)   - tbefore;
            se = TrialStart(correctTrials(tr)) + alignmentSession(correctTrials(tr), end) + tafter;
            trialSpikes(tr) = length(find(spiketrain>sb & spiketrain<se));
        end
        if max(trialSpikes) == 0
            continue   
        end
        
        % Checking stability
        % THIS IS ALL COMMENTED OUT -- NO STABILITY CORRECTION
        baselineSpikes = zeros(1, length(correctTrials));
        for tr = 1:length(correctTrials)
            so = TrialStart(correctTrials(tr)) + alignmentSession(correctTrials(tr), 1);
            baselineSpikes(tr) = length(find(spiketrain>so-0.5 & spiketrain<so));
        end
        B = 10;
        baselineSpikesBatches = reshape(baselineSpikes(1:floor(length(correctTrials)/B)*B), B, []);
        p = anova1(baselineSpikesBatches, [], 'off');
        if p < 0.001
%            display(['Unstable recording detected: ' sessions{s} filenames{f}])
            unstableUnits = unstableUnits + 1;
            
%             tmpFig = figure;
%             hold on
%             boxplot(baselineSpikesBatches)
%             xlabel('Batches of 10 trials')
%             title('Number of spikes in the baseline 0.5 sec')
%             ylabel('Number of spikes')
%             
%             pp = ttest(baselineSpikesBatches);
%             pp(isnan(pp))=0;
%             pp1 = find(diff([0 pp])==1);
%             pp2 = find(diff([pp 0])==-1);
%             [~, pppos] = max(pp2-pp1);
%             ppStart = pp1(pppos);
%             ppEnd = pp2(pppos);
%             
%             if ~isempty(ppStart)
%                 plot(find(pp), min(min(baselineSpikesBatches)), 'b*', 'MarkerSize', 10)
%                 plot(ppStart, min(min(baselineSpikesBatches)), 'r*', 'MarkerSize', 20)
%                 plot(ppEnd, min(min(baselineSpikesBatches)), 'r*', 'MarkerSize', 20)
%             end
%             
%             if isempty(ppStart)
%                 continue
%             end
%             ppStart = (ppStart-1) * B + 1;
%             ppEnd = ppEnd * B;
%             if length(correctTrials) - ppEnd <= B
%                 ppEnd = length(correctTrials);
%             end
%             stableRange = ppStart:ppEnd;
            
            unstableNeuronsMask(unit+1) = 1;
%                       
%             pause
%             close(tmpFig)
%             continue
        else
            %stableRange = 1:length(correctTrials); 
            unstableNeuronsMask(unit+1) = 0;
        end

        % set stable range to the whole session
        % (i.e. no stability correction)
        stableRange = 1:length(correctTrials);
      
        % for debugging purposes: raster plot of this neuron throughout the
        % session
        if 1==0
            tmpFig = figure('Position', [100 700 1600 400]);
            subplot(311)
            hold on
            axis([TrialStart(1) TrialStart(end)+10 -1 2])
            plot(spiketrain, rand(size(spiketrain)), '.')
            for i=1:length(TrialStart)
                plot([TrialStart(i) TrialStart(i)], [-1 2], 'k')
            end
            title([num2str(unit) ': ' sessions{s} filenames{f}], 'interpreter', 'none')
            subplot(312)
            hold on
            plot(TrialStart(correctTrials), stimulus, '*')
            axis([TrialStart(1) TrialStart(end)+10 0 8])
            title('Stimulus id for completed trials only')
            subplot(313)
            h = histc(spiketrain, TrialStart);
            h = h(1:end-1)' ./ diff(TrialStart);
            plot(TrialStart(correctTrials), h(correctTrials), 'r*')
            hold on
            axis([TrialStart(1) TrialStart(end)+10 0 100])
            pause
            close(tmpFig)
        end

        unit = unit+1;
        usedFiles{unit} = [sessions{s} filenames{f}];
        ratMask(unit) = sessionsRats(s);
        pureOdoursMask(unit) = sessionsPure(s);
        
        % Checking for difference between conditions at baseline
%         for st = 1:2
%             for dec = 1:2
%                 trialsSubset = intersect(find(stimulusCategory==st & decision==dec), stableRange);
%                 group(trialsSubset) = (st-1)*2 + dec;
%             end
%         end
%         p = anova1(baselineSpikes(stableRange), group(stableRange), 'off');
%         diffBetweenConditions(unit) = p;
% %          if p < 0.001
% % %             display(['Recording with signif diff at baseline: ' sessions{s} filenames{f}])
% %              anova1(baselineSpikes(stableRange), group(stableRange));
% %              pause
% % %             continue
% %          end
                                  
        if ~isempty(strfind(datadir, 'Feierstein'))
           MaxSt = 2;
           columnSt = stimulusCategory;
        end
        if ~isempty(strfind(datadir, 'Kepecs'))
           MaxSt = 6;
           columnSt = stimulusRatio;
        end

        % loop over possible stimuli / stimulus categories / reward / odour ratios
        for st = 1:MaxSt                              %%%%%
            % loop over possible decisions             %%%
            for dec = 1:2                               % 
                trialsSubset = intersect(find(columnSt==st & decision==dec), stableRange);
                
                % if there is no such combination in this session
                if isempty(trialsSubset)
                    rate(unit, st, dec, :) = nan(1,length(time));
                    rateSTD(unit, st, dec, :) = nan(1,length(time));
                    rateN(unit, st, dec) = 0;    
                    for ii=1:10
                        rateNoise(unit, st, dec, :, ii) = nan(1,length(time));
                    end
                else
                    % loop over correct trials belonging to this
                    % stimulus/decision combination
                    
                    psths = zeros(length(trialsSubset), length(time));
                    
                    for t = 1:length(trialsSubset)
                        trial = correctTrials(trialsSubset(t));
                        
                        tbeg = alignmentSession(trial,1) - tbefore;
                        tend = alignmentSession(trial,end) + tafter;
                        tt = tbeg:0.01:tend;
                        timepoints = [tbeg alignmentSession(trial,:) tend];
                        
                        spiketr = histc(spiketrain - TrialStart(trial), tt);
                        psth = conv(spiketr, gaussKernel, 'same');
                        
                        psthStretched = psths(t,:);
                        
                        % restretching
                        if(restretch)
                            psthStretched = zeros(1,length(time));
                            alignOld = [tbeg alignmentSession(trial,:) tend];
                            alignNew = [-tbefore alignmentMedian alignmentMedian(end)+tafter];

                            for k=1:length(alignNew)-1
                                indtofill = find(time>=alignNew(k) & time<alignNew(k+1));
                                timeint = alignOld(k) + (time(indtofill)-alignNew(k))./ ...
                                    (alignNew(k+1)-alignNew(k))*(alignOld(k+1) - alignOld(k));
                                psthStretched(indtofill) = interp1(tt, psth, timeint);
                            end
                            psthStretched(end) = psth(end);
                        elseif(piecewiseStitch)
                        
                        % piecewise stiching (instead of restretching)
                            psthStretched = [];
                            for event = alignmentSession(trial,:)
                                psthStretched = [psthStretched psth(find(tt>event,1)-45:find(tt>event,1)+45)'];
                            end                        
                        end
                        psths(t,:) = psthStretched;
                    end
                    
                    rate(unit, st, dec, :) = mean(psths, 1);
                    rateSTD(unit, st, dec, :) = std(psths, [], 1);
                    N = size(psths, 1);
                    rateN(unit, st, dec, :) = N;
                    
                    for trtr = 1:size(psths,1)
                        rateAllTrials(unit, st, dec, :, trtr) = psths(trtr,:);
                    end
                end
            end
        end
    end
end        
%%
%display(['Number of unstable recordings found: ' num2str(unstableUnits)])
display(['Number of cells: ' num2str(unit)])

% renaming the variables according to my new convention
timeEvents = alignmentMedian;
timeEventsNames = {'OdorPokeIn','OdorPokeOut','WaterPokeIn','WaterValveOn','WaterPokeOut'};
firingRatesPerTrial_size = size(rateAllTrials);
firingRatesPerTrial_sparse = sparse(rateAllTrials(:));
numOfTrials = rateN;
firingRatesAverage = rate;

%saving
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,2:5,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage(:,2:5,:,:), D, []), 2);
n = find(minN >= 2 & meanFiringRate < 50);
t = 1:length(time);
rateAllTrials= rateAllTrials(n,:,:,:,:);
% save(sprintf('~/gitPublic/elife2016dpca/%s',outputFileName), 'firingRatesPerTrial_sparse', 'firingRatesPerTrial_size', ...
%                      'firingRatesAverage', 'numOfTrials', ...
%                      'time', 'timeEvents', 'timeEventsNames', ...
%                      'unstableNeuronsMask', 'ratMask', 'pureOdoursMask', 'usedFiles','restretch','piecewiseStitch','rateAllTrials');
                 

%%   KENNETH's MAIN PROCESSSING BIT

% load('~/gitPublic/elife2016dpca/data_adam.mat')
D = size(numOfTrials,1);
minN = min(reshape(numOfTrials(:,2:5,:), D, []), [], 2);
meanFiringRate = mean(reshape(firingRatesAverage(:,2:5,:,:), D, []), 2);
n = find(minN >= 2 & meanFiringRate < 50);
t = 1:length(time);

X_all = firingRatesAverage(n,2:5,:,t);

% baseline = mean(mean(mean(X_all,2),3),4);
baseline = nanmean(nanmean(nanmean(firingRatesAverage(n,:,:,:),2),3),4);
X_all = X_all-baseline;

C = size(X_all,2);
D = size(X_all,3);
T = size(X_all,4);
N = size(X_all,1);
fprintf('Loading %d cells.\n',N);

X       = zeros( C*D*T,N);
X_gamma_0 = zeros( C*D*T,N,7);


X_star    = zeros( 2*T,N);
X_star((1:T)+0*T,:) = squeeze(firingRatesAverage(n,1,2,t))';
X_star((1:T)+1*T,:) = squeeze(firingRatesAverage(n,6,1,t))';
X_star = X_star-baseline';

X_t = mean(mean(X_all,2),3);
X_s = mean(mean(X_all,3),4);
X_d = mean(mean(X_all,2),4);
X_st = mean(X_all - X_t - X_s - X_d,3);
X_dt = mean(X_all - X_t - X_s - X_d,2);
X_sd = mean(X_all - X_t - X_s - X_d,4);

X_sdt = ((((((X_all - X_t) - X_s) - X_d) - X_st) - X_dt) - X_sd);



idxs = cell(C,D);

for ii = 1:C
    for jj = 1:D
        idxs{ii,jj} =  (ii-1)*(D*T) + (jj-1)*(T) + (1:T);
        for kk = 1:T
            rr = (ii-1)*(D*T) + (jj-1)*(T) + kk;
            
            X(rr,:) = X_all(:,ii,jj,kk)';
            
            X_gamma_0(rr,:,1) = X_t(:,1,1,kk);
            X_gamma_0(rr,:,2) = X_s(:,ii,1,1);
            X_gamma_0(rr,:,3) = X_d(:,1,jj,1);
            
            X_gamma_0(rr,:,4) = X_st(:,ii,1,kk);
            X_gamma_0(rr,:,5) = X_dt(:,1,jj,kk);
            X_gamma_0(rr,:,6) = X_sd(:,ii,jj,1);
            
            X_gamma_0(rr,:,7) = X_sdt(:,ii,jj,kk);
            
        end
    end
end
X_gamma = zeros(C*D*T,N,4);
X_gamma(:,:,1) = X_gamma_0(:,:,1);
X_gamma(:,:,2) = sum(X_gamma_0(:,:,[2 4]),3);
X_gamma(:,:,3) = sum(X_gamma_0(:,:,[3 5]),3);
X_gamma(:,:,4) = sum(X_gamma_0(:,:,[6 7]),3);

save('Data/kepecsProcessed.mat','-v7.3','X_gamma','X_gamma_0','X','X_star','T','D','C','idxs','X_all');
            