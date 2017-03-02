%% stimulationExperiment
close all; clear all; clc

%% GENERATE PARAMS
[Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams();
%% change default parameters
T = 5*1000; % shortened trials
note = ''; % if nothing is different
% note = 'break E-PV'; % 

% warning('breaking PV connections!');
% gSyn(1,2) = 0;
%% genweights: generate recurrent and external weights and external input
%  rates
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

%% set up stimulation

% other parameters
NT = round(T/dt);
nTrials = 25;

stimPop = 2; % PV Neurons
curstim = 200; % picoAmps
stimLength = 100; % ms
stimStart = 2.5*1000; % ms

% set up groups
nStimGroups = 2;
groupToStim = 1;
groupInds = pinds(stimPop):pinds(stimPop+1)-1;

% big matrix
stimI = zeros(Ntot,NT);

groupIndsRand = groupInds(randperm(length(groupInds)));

nPerGroup = length(groupInds)/nStimGroups;
for i=1:nStimGroups
    startInd = (i-1)*nPerGroup+1;
    endInd = i*nPerGroup;
    thisGroupInds = groupIndsRand(startInd:endInd);
    
    if i==groupToStim
        stimulatedNeurons = thisGroupInds;
        stimI(thisGroupInds,stimStart/dt:stimStart/dt+stimLength/dt) = ...
            curstim;
    end
    
end
%% initialize other parameters

% state vectors
v = -60.*ones(Ntot,1);
vT = zeros(Ntot,1);
lastSpike = -100.*ones(Ntot,1);
whichpop = zeros(Ntot,1);

for pp=1:Npop
    % makes whichpop a vector where each neuron gets a population id
    whichpop(pinds(pp):pinds(pp+1)-1) = pp;
    % calculate gaussian-random threshold
    vT(pinds(pp):pinds(pp+1)-1) = vTpop(pp) + sigvT(pp)*randn(Ncells(pp),1);
end
nextext = zeros(Ntot,1);
% set time to next external stimulation. Note that with higher external
% firing rates, since rate_ext (rext) is in the denominator, time to next
% one decreases as firing rate increases. exprnd(1) is a random draw from
% the exponential distribution with mean (scale factor) == 1
% ~~~NOTE!!!~~~ I think this translates average input rate to KILOHERTZ
% when dT is set to 0.1 ms! Should rewrite this to take dT into account
for cc = 1:Ntot
    nextext(cc) = exprnd(1)/rext(cc);
end

w_adapt = zeros(1,Ntot);
% for each cell, one set of ODEs for each presynaptic population
xrise = zeros(4,Ntot);
xdecay = zeros(4,Ntot);
% for each cell, one set of ODEs for each postsynaptic population
D = ones(4,Ntot);
F = ones(4,Ntot);

% return values, not sure if necessary
ns = 0;
% changed this to not be .05 times T
maxns = round(T*Ntot);
% maxns = inf;
times = zeros(1,maxns);
tinds = zeros(1,maxns);
vall = zeros(1);

%% JIM ADDED
allT = (dt:dt:T);
allV = nan(Ntot,NT);
allF = nan(Ntot,NT);
allD = nan(Ntot,NT);
allRise = nan(Ntot,NT);
allDecay = nan(Ntot,NT);
allExt = zeros(Ntot,NT);
allSpikes = zeros(Ntot,NT);

% 100 = 10Hz imaging, 40 = 25 Hz
binSize = 40; % in milliseconds

[~,downsampledT] = downsampleSpikes(allSpikes(1,:),binSize,dt);

fRates = nan(Ntot,length(downsampledT),nTrials);
%% SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('starting sim');
for trial = 1:nTrials
    % generate new weight matrices
    [rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);
    
    % make sure important state vectors are reset
    v = -60.*ones(Ntot,1);
    lastSpike = -100.*ones(Ntot,1);
    nextext = zeros(Ntot,1);
    % set time to next external stimulation. Note that with higher external
    % firing rates, since rate_ext (rext) is in the denominator, time to next
    % one decreases as firing rate increases. exprnd(1) is a random draw from
    % the exponential distribution with mean (scale factor) == 1
    % ~~~NOTE!!!~~~ I think this translates average input rate to KILOHERTZ
    % when dT is set to 0.1 ms! Should rewrite this to take dT into account
    for cc = 1:Ntot
        nextext(cc) = exprnd(1)/rext(cc);
    end
    w_adapt = zeros(1,Ntot);
    % for each cell, one set of ODEs for each presynaptic population
    xrise = zeros(4,Ntot);
    xdecay = zeros(4,Ntot);
    % for each cell, one set of ODEs for each postsynaptic population
    D = ones(4,Ntot);
    F = ones(4,Ntot);
    % changed this to not be .05 times T
    maxns = round(T*Ntot);
    % maxns = inf;
    times = zeros(1,maxns);
    tinds = zeros(1,maxns);
    vall = zeros(1);   
    
    for tt = 1:NT
        t = dt*tt;
        
        if mod(tt,NT/100)==1 % print percent complete
            fprintf('%d%%\n', round(100*tt/NT));
        end
        
        % update synaptic / adaptation parameters
        for cc=1:Ntot
            pc = whichpop(cc);
            % update w term (adaptve part of aEIF model, brette and gerstner)
            w_adapt(cc) = w_adapt(cc) + ...
                (dt/tauw_adapt(pc))*(a_adapt(pc)*(v(cc)-EL(pc))-w_adapt(cc));
            if isnan(w_adapt(cc))
                disp('w error');
            end
            for qq=1:Npop
                xrise(qq,cc) = xrise(qq,cc)- dt*xrise(qq,cc)/taurise(qq,pc);
                xdecay(qq,cc) = xdecay(qq,cc) - dt*xdecay(qq,cc)/taudecay(qq,pc);
                D(qq,cc) = D(qq,cc) + dt*(1-D(qq,cc))/tauD(pc,qq);
                F(qq,cc) = F(qq,cc) + dt*(1-F(qq,cc))/tauF(pc,qq);
            end
        end
        
        for cc = 1:Ntot
            pc = whichpop(cc);
            
            if t> (lastSpike(cc)+tauref(pc)) % not in refractory period
                % aEIF MODEL IMPLEMENTATION HERE!!!
                % external current - adaptation current - leak current ...
                % + exponential current for spike initiation
                dv = stimI(cc,tt) - w_adapt(cc)-gL(pc)*(v(cc)-EL(pc))+...
                    gL(pc)*deltaT(pc)*exp((v(cc)-vT(cc))/deltaT(pc));
                
                % ADD SYNAPTIC CURRENTS HERE! Current depends on population
                for qq=1:Npop
                    % minus sign because if xDecay < xRise, term is negative
                    dv = dv - gSyn(qq,pc)*(v(cc)-Esyn(qq,pc))*(xdecay(qq,cc)-xrise(qq,cc))/(taudecay(qq,pc)-taurise(qq,pc));
                end
                
                %             if dv == inf
                %                 disp('dv2 error');
                %             end
                if dv == inf
                    v(cc) = 20;
                else
                    % change current to voltage
                    v(cc) = v(cc)+dt*dv/Cm(pc);
                end
                if isnan(v(cc))
                    disp('error');
                end
            end % end if in refractory period
        end
        
        % DO SPIKES!
        for cc = 1:Ntot
            pc = whichpop(cc);
            
            % not sure why max spike ns is even a thing
            if (v(cc)>vth(pc))&&(ns<maxns) % spike occurred
                v(cc) = vre(pc);
                lastSpike(cc) = t;
                w_adapt(cc) = w_adapt(cc) + b_adapt(pc);
                ns = ns+1;
                times(ns) = t;
                tinds(ns) = cc;
                
                % propagate spike
                % wind contains ids of downstream targets
                for kk=wind(cc):wind(cc+1)-1
                    ipost = wipost(kk); % index of downstream target
                    ppost = whichpop(ipost);
                    xrise(pc,ipost) = xrise(pc,ipost) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                    xdecay(pc,ipost) = xdecay(pc,ipost) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                end
                
                for qq = 1:Npop
                    F(qq,cc) = F(qq,cc)+ UF(pc,qq)*(Fmax(pc,qq)-F(qq,cc));
                    D(qq,cc) = D(qq,cc)*UD(pc,qq);
                end
                
                allSpikes(cc,tt) = 1;
            end
        end
        
        % do external excitation. NOTE as above!! This give kHz external
        % stimulation, doesn't vary with dT!!!
        for cc=1:Ntot
            while(t>nextext(cc))
                nextext(cc) = nextext(cc)+ exprnd(1)/rext(cc);
                % does this mean external excitation only goes to E neurons?
                xrise(1,cc) = xrise(1,cc)+wext(cc); % increment excitation
                xdecay(1,cc) = xdecay(1,cc)+wext(cc);
                allExt(cc,tt) = allExt(cc,tt)+1;
            end
        end
        
        % JIM ADDED
        allV(:,tt) = v;
        allF(:,tt) = mean(F)'; % average over all input pops
        allD(:,tt) = mean(D)';
        allRise(:,tt) = mean(xrise)';
        allDecay(:,tt) = mean(xdecay)';
    end % end sim
    
    disp('100%');
    disp('Calculating rates');
    counts = zeros(1,Ntot);
    for cc=1:Ntot
        counts(cc) = sum(tinds==cc);
    end
    
    rates = 1000.*counts/(T);
    
    for pp=1:Npop
        fprintf('Population %d firing rate: %.2f\n', pp,...
            mean(rates(pinds(pp):pinds(pp+1)-1)));
    end
    
    times = times(1:ns);
    tinds = tinds(1:ns);
    
    
    disp('Binning firing rates');
    for cc=1:Ntot
        [fRate,~] = downsampleSpikes(allSpikes(cc,:),binSize,dt);
        fRates(cc,:,trial) = fRate;
    end
    
    
end
%% raster figure

% figure;
% subplot(4,3,1:9);
% for i=1:Ntot
%     spikes = allSpikes(i,:);
%     spikeTimes = allT(find(spikes));
%     if ~isempty(spikeTimes)
%         if sum(stimulatedNeurons==i)
%             for j=1:length(spikeTimes)
%                 line([spikeTimes(j) spikeTimes(j)], [i-1 i], 'color', ...
%                     'b');
%             end
%         else
%             for j=1:length(spikeTimes)
%                 line([spikeTimes(j) spikeTimes(j)], [i-1 i], 'color', ...
%                     'k');
%             end
%         end
%     end
% end
% ylim([0 Ntot]);
% subplot(4,3,10:12);
% [row,col] = find(stimI==max(max(stimI,[],2)),1);
% plot(allT,stimI(row,:));
% 
% %% average firing rate figure
% binSize = 100; % in milliseconds
% groupNames = {'E', 'PV', 'SOM', 'VIP', 'stimE'};
% [~,downsampledT] = downsampleSpikes(allSpikes(1,:),binSize,dt);
% 
% fRates = nan(Ntot,length(downsampledT));
% disp('Binning firing rates');
% for cc=1:Ntot
%     [fRate,~] = downsampleSpikes(allSpikes(cc,:),binSize,dt);
%     fRates(cc,:) = fRate;
% end
% 
% groupAverages = zeros(nStimGroups+1,length(downsampledT));
% gNums = zeros(nStimGroups+1,1);
% gNum = 1;
% for i=1:Ntot
%     if sum(stimulatedNeurons==i)>0
%         groupAverages(nStimGroups+1,:) = groupAverages(nStimGroups+1,:)+...
%             fRates(i,:);
%         gNums(nStimGroups+1) = gNums(nStimGroups+1)+1;
%     else
%         groupAverages(whichpop(i),:) = groupAverages(whichpop(i),:)+...
%             fRates(i,:);
%         gNums(whichpop(i)) = gNums(whichpop(i))+1;
%     end
% end
% % firing rate averages, in hertz
% groupAverages = groupAverages./gNums.*binSize*dt;
% 
% 
% figure;
% subplot(4,3,1:9);
% hold on;
% for i=1:nStimGroups+1
%     plot(downsampledT,groupAverages(i,:));
% end
% legend(groupNames);
%% firing rate figure for single trial
fR = fRates(:,:,randi(size(fRates,3)));
cMap = [
    0 0.3 0.6; % blue-ish
    0 .5 0; % green
    1 .5 0; % orange
    115/255 44/255 123/255 % purpley
    0 0 0];
% the color of the stimulated population (#5) should be a darker version of
% the normal color
cMap(5,:) = cMap(stimPop,:)./2; 

groupNames = {'E', 'PV', 'SOM', 'VIP'};
groupNames{5} = ['stim',groupNames{stimPop}];
nG = length(groupNames);

% COME UP WITH DESCRIPTIVE FIGURE NAME FOT THIS FIGURE
experimentDescription = sprintf(['StimPop_%s_groupsInPop_%d_photoCurrent_%dpA_'...
    'stimLength_%dms%s'], groupNames{stimPop}, nStimGroups, curstim, stimLength,...
    note);
figName = [experimentDescription '_singleTrialPopFRates'];
saveLoc = 'D:\Analysis\cortexModel\';

groupAverages = zeros(nG,length(downsampledT));
gNums = zeros(nG,1);
for i=1:Ntot
    if sum(stimulatedNeurons==i)>0
        groupAverages(5,:) = groupAverages(5,:)+...
            fR(i,:);
        gNums(5) = gNums(5)+1;
    else
        groupAverages(whichpop(i),:) = groupAverages(whichpop(i),:)+...
            fR(i,:);
        gNums(whichpop(i)) = gNums(whichpop(i))+1;
    end
end
% firing rate averages, in hertz
groupAverages = groupAverages./gNums.*binSize*dt;

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
for i=1:nG
    plotH = plot(downsampledT./1000,groupAverages(i,:));
    plotH.Color = cMap(i,:);
    plotH.LineWidth = 1;
end
% Use Tex to change the color of the text on the figure
legendCell = cell(nG,1);
for i=1:nG
    legendCell{i} = sprintf('\\color[rgb]{%.2f, %.2f, %.2f} %s',...
        cMap(i,:), groupNames{i});
end
legH = legend(legendCell);
legH.FontSize = 16;
% Make figure pretty
a = gca;
a.FontSize = 14;
a.XTick = [];
a.TickDir = 'out';
title(strrep(figName,'_', ' '), 'fontsize', 18);
ylabel('FR (Hz)', 'fontsize', 16);
legH.FontWeight = 'bold';

% Plot stimulation current!
subplot(4,3,10:12);
% find a good zoom window around the stimulation
if curstim > 0
    [row,col] = find(stimI==max(max(stimI,[],2)),1);
elseif curstim < 0
    [row,col] = find(stimI==min(min(stimI,[],2)),1);
end
plotH = plot(allT./1000,stimI(row,:));
plotH.LineWidth = 1;
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);
xlabel('Time (s)', 'fontsize', 16);
export_fig([saveLoc, figName], '-png', '-eps');

%% zoom in figure around stimulation
figName = [experimentDescription '_singleTrialPopFRatesZoomed' note];

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
for i=1:nG
    plotH = plot(downsampledT./1000,groupAverages(i,:));
    plotH.Color = cMap(i,:);
    plotH.LineWidth = 1;
end
% Use Tex to change the color of the text on the figure
legendCell = cell(nG,1);
for i=1:nG
    legendCell{i} = sprintf('\\color[rgb]{%.2f, %.2f, %.2f} %s',...
        cMap(i,:), groupNames{i});
end
legH = legend(legendCell);
legH.FontSize = 16;
% Make figure pretty
a = gca;
a.FontSize = 14;
a.XTick = [];
a.TickDir = 'out';
title(strrep(figName,'_', ' '), 'fontsize', 18);
ylabel('FR (Hz)', 'fontsize', 16);
legH.FontWeight = 'bold';

% find a good zoom window around the stimulation
if curstim > 0
    [row,col] = find(stimI==max(max(stimI,[],2)),1);
    dI = diff(stimI(row,:));
    [~,stimOn ] = max(dI);
    [~,stimOff ] = min(dI);
    [val,downsampledStartInd] = min(abs(downsampledT-allT(stimOn)));
    [val,downsampledEndInd] = min(abs(downsampledT-allT(stimOff)));
elseif curstim < 0
    [row,col] = find(stimI==min(min(stimI,[],2)),1);
    dI = diff(stimI(row,:));
    [~,stimOn ] = min(dI);
    [~,stimOff ] = max(dI);
    [val,downsampledStartInd] = min(abs(downsampledT-allT(stimOn)));
    [val,downsampledEndInd] = min(abs(downsampledT-allT(stimOff)));
end

zoomWidth = 15; % in bins
% ZOOM IN! 
xlim([downsampledT(downsampledStartInd-zoomWidth)./1000 ...
    downsampledT(downsampledStartInd+zoomWidth)./1000]);

% Plot stimulation current!
subplot(4,3,10:12);
downsampledStim = downsample(stimI(row,:),binSize/dt);
plotH = plot(downsampledT./1000,downsampledStim);
plotH.LineWidth = 1;
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);
xlim([downsampledT(downsampledStartInd-zoomWidth)./1000 ...
    downsampledT(downsampledStartInd+zoomWidth)./1000]);
xlabel('Time (s)', 'fontsize',16);

export_fig([saveLoc, figName], '-png', '-eps');

%% AVERAGE ALL TRIALS TOGETHER
figName = sprintf('%s_%dTrialsAvgPopRates%s', experimentDescription, nTrials,note);

fR = mean(fRates,3);

groupAverages = zeros(nG,length(downsampledT));
gNums = zeros(nG,1);
gNum = 1;
for i=1:Ntot
    if sum(stimulatedNeurons==i)>0
        groupAverages(nG,:) = groupAverages(nG,:)+...
            fR(i,:);
        gNums(nG) = gNums(nG)+1;
    else
        groupAverages(whichpop(i),:) = groupAverages(whichpop(i),:)+...
            fR(i,:);
        gNums(whichpop(i)) = gNums(whichpop(i))+1;
    end
end
% firing rate averages, in hertz
groupAverages = groupAverages./gNums.*binSize*dt;

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
for i=1:nG
    plotH = plot(downsampledT./1000,groupAverages(i,:));
    plotH.Color = cMap(i,:);
    plotH.LineWidth = 1;
end
% Use Tex to change the color of the text on the figure
legendCell = cell(nG,1);
for i=1:nG
    legendCell{i} = sprintf('\\color[rgb]{%.2f, %.2f, %.2f} %s',...
        cMap(i,:), groupNames{i});
end
legH = legend(legendCell);
legH.FontSize = 16;
% Make figure pretty
a = gca;
a.FontSize = 14;
a.XTick = [];
a.TickDir = 'out';
title(strrep(figName,'_', ' '), 'fontsize', 18);
ylabel('FR (Hz)', 'fontsize', 16);
legH.FontWeight = 'bold';

% Plot stimulation current!
subplot(4,3,10:12);
% [row,col] = find(stimI==max(max(stimI,[],2)),1);
plot(allT./1000,stimI(row,:));
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);

export_fig([saveLoc, figName], '-png', '-eps');
%% zoom in figure around stimulation
figName = sprintf('%s_%dTrialsAvgPopRatesZoomed%s', ...
    experimentDescription, nTrials,note);

fR = mean(fRates,3);

groupAverages = zeros(nG,length(downsampledT));
gNums = zeros(nG,1);
gNum = 1;
for i=1:Ntot
    if sum(stimulatedNeurons==i)>0
        groupAverages(5,:) = groupAverages(5,:)+...
            fR(i,:);
        gNums(5) = gNums(5)+1;
    else
        groupAverages(whichpop(i),:) = groupAverages(whichpop(i),:)+...
            fR(i,:);
        gNums(whichpop(i)) = gNums(whichpop(i))+1;
    end
end
% firing rate averages, in hertz
groupAverages = groupAverages./gNums.*binSize*dt;

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
for i=1:nG
    plotH = plot(downsampledT./1000,groupAverages(i,:));
    plotH.Color = cMap(i,:);
    plotH.LineWidth = 1;
end
% Use Tex to change the color of the text on the figure
legendCell = cell(nG,1);
for i=1:nG
    legendCell{i} = sprintf('\\color[rgb]{%.2f, %.2f, %.2f} %s',...
        cMap(i,:), groupNames{i});
end
legH = legend(legendCell);
legH.FontSize = 16;
% Make figure pretty
a = gca;
a.FontSize = 14;
a.XTick = [];
a.TickDir = 'out';
title(strrep(figName,'_', ' '), 'fontsize', 18);
ylabel('FR (Hz)', 'fontsize', 16);
legH.FontWeight = 'bold';
xlim([downsampledT(downsampledStartInd-zoomWidth)/1000 ...
    downsampledT(downsampledStartInd+zoomWidth)/1000]);

% Plot stimulation current!
subplot(4,3,10:12);
plot(allT./1000,stimI(row,:));
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);
xlim([downsampledT(downsampledStartInd-zoomWidth)/1000 ...
    downsampledT(downsampledStartInd+zoomWidth)/1000]);

export_fig([saveLoc, figName], '-png', '-eps');

%% Average of fR 1 second before, compare to absolute peak fR during stim
figName = sprintf('%s_%dTrials_AvgRateBeforeDuringStim%s', ...
    experimentDescription, nTrials,note);

avgBefore = cell(nG,1);
avgDuring = cell(nG,1);
for cc = 1:Ntot
    
    % get the avg firing rate across all trials 1 second before stim
    fRateBefore = fR(cc,downsampledStartInd-1000/binSize:downsampledStartInd-1);
    % get avg firing rate during all trials during stim
    fRateDuring = fR(cc,downsampledStartInd:downsampledEndInd);
    
    % if stimulated population
    if sum(cc==stimulatedNeurons)>0
        % add to cell array
        avgArray = avgBefore{5};
        avgArray(end+1) = mean(fRateBefore);
        avgBefore{5} = avgArray;
        
        avgArray = avgDuring{5};
        avgArray(end+1) = mean(fRateDuring);
        avgDuring{5} = avgArray;
    else
        % add to cell array of specified population
        avgArray = avgBefore{whichpop(cc)};
        avgArray(end+1) = mean(fRateBefore);
        avgBefore{whichpop(cc)} = avgArray;
        
        avgArray = avgDuring{whichpop(cc)};
        avgArray(end+1) = mean(fRateDuring);
        avgDuring{whichpop(cc)} = avgArray;
    end
end

meanBefore = cellfun(@mean,avgBefore);
meanDuring = cellfun(@mean,avgDuring);

% make new cell array with distributions side-by-side
frSideBySide = cell(nG*2,1);
inds = 1:2:nG*2;
% make a new colormap to match this array
cMap2 = nan(nG*2,3);
for i=1:nG
    frSideBySide{inds(i)} = avgBefore{i};
    frSideBySide{inds(i)+1} = avgDuring{i};
    cMap2(inds(i),:) = cMap(i,:);
    cMap2(inds(i)+1,:) = cMap(i,:);
end

meanSideBySide = cellfun(@mean,frSideBySide);

% figure out the maximum firing rate to set the y limit
allFR = [avgBefore{:} avgDuring{:}];
maxFR = ceil(max(allFR));

figure('Units', 'Pix', 'outerposition', [0 50 1600 1150], 'color', [1 1 1]);
plotSpread(frSideBySide, 'distributionColors', cMap2);
hold on;
for i=1:nG*2
    line([i-.4 i+.4], [meanSideBySide(i) meanSideBySide(i)], 'Color', ...
        cMap2(i,:), 'lineWidth', 1.5);
end
hold off;
newLegendCell = cell(nG*2,1);
for i=1:nG
    newLegendCell{inds(i)} = legendCell{i};
    newLegendCell{inds(i)+1} = '';
end
legH = legend(newLegendCell);
legH.FontSize = 16;
a = gca;
a.YLim = [0 maxFR];
a.FontSize = 14;
ylabel('Firing rate (Hz)', 'fontsize', 18);
a.XTick = [];

title(strrep(figName,'_', ' '), 'Fontsize', 18,'fontweight', 'bold');
export_fig([saveLoc, figName], '-png', '-eps');