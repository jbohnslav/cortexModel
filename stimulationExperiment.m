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
% note = 'more SOM bg'; % 

% warning('breaking PV connections!');
% gSyn(1,2) = 0;
% r0(1) = r0(1)*2;
% r0(3) = r0(3)*2;
% gSyn(1,1) = gSyn(1,1)*2;
% r0(3) = r0(3)*2;
% Ncells = Ncells.*4;
%% genweights: generate recurrent and external weights and external input
%  rates
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

%% set up stimulation

% other parameters
NT = round(T/dt);
nTrials = 1;
onePulse = 1;

stimPop = 4; % 
curstim = 250; % picoAmps
stimLength = 50; % ms
stimStart = 2.5*1000; % ms

% set up groups
nStimGroups = 4;
groupToStim = 1;
groupInds = pinds(stimPop):pinds(stimPop+1)-1;

% big matrix
stimI = zeros(Ntot,NT);

groupIndsRand = groupInds(randperm(length(groupInds)));
nPerGroup = length(groupInds)/nStimGroups;
stimWave = zeros(1,NT);
allT = (dt:dt:T);
if onePulse
    stimTInds = stimStart/dt:stimStart/dt+stimLength/dt-1;
    stimWave(stimTInds) = curstim;
else
    % multiple pulses at a given pulse rate
    Hz = 20;
    stimLength = 50; % ms
    nPulses = 10;
    totalStimLengthT = nPulses/Hz; % seconds
    stimTInds = stimStart/dt:stimStart/dt+totalStimLengthT*1000/dt-1;
    % duty cycle is the percent of the period where the square wave is +
    % 50 is standard
    dutyCycle = 50;
    squareWave = square(allT(stimTInds).*(2*pi/1000*Hz),dutyCycle);
    squareWave(squareWave<0) = 0; % rectify! not below zero
    stimWave(stimTInds) = squareWave.*curstim;
end

for i=1:nStimGroups
    startInd = round((i-1)*nPerGroup+1);
    endInd = round(i*nPerGroup);
    thisGroupInds = groupIndsRand(startInd:endInd);
    
    if i==groupToStim
        stimulatedNeurons = thisGroupInds;
        stimI(stimulatedNeurons,:) = repmat(stimWave,[length(thisGroupInds) 1]);
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

allV = nan(Ntot,NT,nTrials);
allF = nan(Ntot,NT);
allD = nan(Ntot,NT);
allRise = nan(Ntot,NT);
allDecay = nan(Ntot,NT);
allExt = zeros(Ntot,NT);
allSpikes = zeros(Ntot,NT,nTrials);

% 100 = 10Hz imaging, 40 = 25 Hz
binSize = 40; % in milliseconds

% [~,downsampledT] = downsampleSpikes(allSpikes(1,:,1),binSize,dt);

fRates = nan(Ntot,NT,nTrials);
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
                
                allSpikes(cc,tt,trial) = 1;
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
        allV(:,tt,trial) = v;
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
%         [fRate,~] = downsampleSpikes(allSpikes(cc,:,trial),binSize,dt);
%         gauss_window = 100/dt; % 100 ms bin width for downsampling
        gauss_SD = 25; % ms standard deviation
        edges = -3*gauss_SD:dt:3*gauss_SD;
        kernel = normpdf(edges,0,gauss_SD);
        sdf = conv(allSpikes(cc,:,trial),kernel,'same').*1000;
        fRates(cc,:,trial) = sdf;
    end
    
    
end
%% raster figure
disp('making raster fig...');
cMap = [
    0 0.3 0.6; % blue-ish
    0 .5 0; % green
    1 .5 0; % orange
    115/255 44/255 123/255 % purpley
    0 0 0];
cMap(5,:) = cMap(stimPop,:)./2;
groupNames = {'E', 'PV', 'SOM', 'VIP'};
groupNames{5} = ['stim',groupNames{stimPop}];
nG = length(groupNames);
% 
% % COME UP WITH DESCRIPTIVE FIGURE NAME FOT THIS FIGURE
experimentDescription = sprintf(['StimPop_%s_nGroups_%d_photoCurrent_%dpA_'...
    'stimLength_%dms_%s'], groupNames{stimPop}, nStimGroups, curstim, stimLength,...
    note);
figName = [experimentDescription '_singleTrialRaster'];
saveLoc = 'D:\Analysis\cortexModel\';

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
trial = randi(nTrials);
for i=1:Ntot
    spikes = allSpikes(i,:,trial);
    spikeTimes = allT(find(spikes));
    if ~isempty(spikeTimes)
        if sum(stimulatedNeurons==i)
               plot(spikeTimes,ones(length(spikeTimes),1).*i,'.',...
                   'Color', cMap(5,:));
        else
            for j=1:length(spikeTimes)
                plot(spikeTimes,ones(length(spikeTimes),1).*i,'.',...
                   'Color', cMap(whichpop(i),:));
            end
        end
    end
    fprintf('cc: %d\n', i);
end
hold off;
ylim([0 Ntot]);

legendCell = cell(nG,1);
for i=1:nG
    legendCell{i} = sprintf('\\color[rgb]{%.2f, %.2f, %.2f} %s',...
        cMap(i,:), groupNames{i});
end
legH = legend(legendCell);
legH.FontSize = 16;
legH.FontWeight = 'bold';

a = gca;
a.FontSize = 14;
a.XTick = [];
a.TickDir = 'out';
title(strrep(figName,'_', ' '), 'fontsize', 18);
ylabel('Cell', 'fontsize', 16);


subplot(4,3,10:12);
plot(allT./1000,stimWave);
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);
xlabel('Time (s)', 'fontsize', 16);
disp('exporting figure...');
export_fig([saveLoc, figName], '-png');

%% AVERAGE ALL TRIALS TOGETHER
figName = sprintf('%s_%dTrialsAvgPopRates%s', experimentDescription, nTrials);

fR = mean(fRates,3);

groupAverages = zeros(nG,length(allT));
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
groupAverages = groupAverages./gNums;

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
for i=1:nG
    plotH = plot(allT./1000,groupAverages(i,:));
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
plot(allT./1000,stimWave);
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);

export_fig([saveLoc, figName], '-png');
%% zoom in figure around stimulation
zoomWidth = 100; % milliseconds
stimOn = allT(stimTInds(1));
stimOff = allT(stimTInds(end)+1);

figName = sprintf('%s_%dTrialsAvgPopRatesZoomed%s', ...
    experimentDescription, nTrials);

% work with means for each neuron across trials
fR = mean(fRates,3); 

% make a cell array of matrices, each matrix has fR x cells
groupAverages = cell(nG,1);
nInGroup = nan(nG,1);
for i=1:length(Ncells)
    if i==stimPop
        nInGroup(i) = Ncells(i)- length(stimulatedNeurons);
        nInGroup(5) = length(stimulatedNeurons);
    else
        nInGroup(i) = Ncells(i);
    end
end
% pre-allocate for speed
for i=1:nG
    groupAverages{i} = nan(nInGroup(i),NT);
end
howMany = ones(nG,1);
for i=1:Ntot
    % if this neuron is in stimulated group
    if sum(stimulatedNeurons==i)>0
%         mat = groupAverages{5};
        groupAverages{5}(howMany(5),:) = fR(i,:);
        howMany(5) = howMany(5)+1;
    else
%         mat = groupAverages{whichpop(i)};
        groupAverages{whichpop(i)}(howMany(whichpop(i)),:) = fR(i,:);
        howMany(whichpop(i)) = howMany(whichpop(i))+1;
    end
%     i
end
% firing rate averages, in hertz
groupAvgAllCells = nan(nG,NT);
groupStdAllCells = nan(nG,NT);
for i=1:nG
    groupAvgAllCells(i,:) = mean(groupAverages{i});
    groupStdAllCells(i,:) = std(groupAverages{i});
end


figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);
subplot(4,3,1:9);
hold on;
for i=1:nG
    SEB = shadedErrorBar(allT./1000,groupAvgAllCells(i,:),...
        groupStdAllCells(i,:)/sqrt(nInGroup(i)));
    SEB.mainLine.Color = cMap(i,:);
    SEB.mainLine.LineWidth = 1.5;
    SEB.patch.FaceAlpha = .3;
    SEB.patch.FaceColor = cMap(i,:);
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
xlim([(stimOn-zoomWidth)/1000 (stimOff+zoomWidth)/1000]);

% Plot stimulation current!
subplot(4,3,10:12);
plot(allT./1000,stimWave);
a = gca;
a.FontSize = 14;
a.TickDir = 'out';
ys = ylim;
ylim([ys(1), ys(2)+50]);
ylabel({'External stimulation', 'current (pA)'}, 'fontsize', 16);
xlim([(stimOn-zoomWidth)/1000 (stimOff+zoomWidth)/1000]);

export_fig([saveLoc, figName], '-png', '-eps');

%% Average of fR 1 second before, compare to absolute peak fR during stim
figName = sprintf('%s_%dTrials_AvgRateBeforeDuringStim%s', ...
    experimentDescription, nTrials);

avgBefore = cell(nG,1);
avgDuring = cell(nG,1);
for cc = 1:Ntot
    
    % get the avg firing rate across all trials 100 ms before stim
    fRateBefore = fR(cc,stimTInds(1)-100/dt:stimTInds(1)-1);
    % get avg firing rate during / after all trials during stim
    fRateDuring = fR(cc,stimTInds(1):stimTInds(1)+100/dt);
    
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
export_fig([saveLoc, figName], '-png');

%%
figName = [experimentDescription, 'frRateCDF'];

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'color', [1 1 1]);

subplot(1,2,1);
scH = scatter(frSideBySide{1},frSideBySide{2}, 30,'filled');
scH.MarkerFaceColor = cMap(1,:);
hold on;
% axis([0 1 0 1])
maxFR = max([frSideBySide{1} frSideBySide{2}]);
axis([0 maxFR 0 maxFR]);
axis square;
xs = xlim;
ys = ylim;
line(xs,ys,'color', 'k');
mdl = fitlm(frSideBySide{1},frSideBySide{2});
coeffs = mdl.Coefficients.Estimate;
coeffs(mdl.Coefficients.pValue>.05) = 0;
x = linspace(xs(1),xs(2),1000);
y = x*coeffs(2)+coeffs(1);
plotH = plot(x,y);

axis square;
grid on;
a = gca;
a.FontSize = 14;
legend('Data','y=x', sprintf('y=%.2fx + %.2f,r^2=%.3f', coeffs(2),coeffs(1),...
    mdl.Rsquared.Ordinary));
xlabel('FR Rate before stim', 'fontsize', 16, 'fontweight', 'bold');
ylabel('FR Rate after stim', 'fontsize', 16, 'fontweight', 'bold'); 

subplot(1,2,2);
frBins = linspace(0,maxFR,25); % Hz
[f,x] = hist(frSideBySide{1},frBins);
normFPreStim = f/sum(f);

[f,x] = hist(frSideBySide{2},frBins);
normFPostStim = f/sum(f);

plotH = plot(frBins,cumsum(normFPreStim));
plotH.Color = cMap(1,:);
plotH.LineWidth = 2;
hold on;
plotH = plot(frBins,cumsum(normFPostStim));
plotH.Color = cMap(1,:)./2;
plotH.LineWidth = 2;


legH = legend('E before stim','E during stim');
legH.FontSize = 16;
a = gca;
a.FontSize = 14;
axis square;
grid on;

xlabel('Firing rate', 'fontsize', 16, 'fontweight', 'bold');
ylabel('CDF', 'fontsize', 16, 'fontweight', 'bold');

[h,p,k] = kstest2(frSideBySide{1},frSideBySide{2},'Tail','Larger');

title(sprintf('Firing Rate CDF: Stim > Before P=%.4f', p));
export_fig([saveLoc, figName], '-png');
%%
% figName = sprintf('%s_%dTrials_AvgVoltage', ...
%     experimentDescription, nTrials);
% 
% groupAverages = zeros(nG,length(allT));
% gNums = zeros(nG,1);
% gNum = 1;
% for i=1:Ntot
%     if sum(stimulatedNeurons==i)>0
%         groupAverages(nG,:) = groupAverages(nG,:)+...
%             mean(allV(i,:,:),3);
%         gNums(nG) = gNums(nG)+1;
%     else
%         groupAverages(whichpop(i),:) = groupAverages(whichpop(i),:)+...
%             mean(allV(i,:,:),3);
%         gNums(whichpop(i)) = gNums(whichpop(i))+1;
%     end
% end
% % firing rate averages, in hertz
% groupAverages = groupAverages./gNums;
% 
% figure('Units', 'Pix', 'outerposition', [0 50 1600 1150], 'color', [1 1 1]);
% subplot(4,1,1:3);
% hold on;
% for i=1:nG
%     plot(allT,groupAverages(i,:),'Color',cMap(i,:));
% end
% hold off;
% legend(groupNames)
% legendCell = cell(nG,1);
% for i=1:nG
%     legendCell{i} = sprintf('\\color[rgb]{%.2f, %.2f, %.2f} %s',...
%         cMap(i,:), groupNames{i});
% end
% legH = legend(legendCell);
% legH.FontSize = 16;
% % Make figure pretty
% a = gca;
% a.FontSize = 14;
% a.XTick = [];
% a.TickDir = 'out';
% title(strrep(figName,'_', ' '), 'fontsize', 18);
% ylabel('V (mV)', 'fontsize', 16);
% legH.FontWeight = 'bold';
% width = 100; %ms
% xlim([allT(stimTInds(1))-width allT(stimTInds(end))+width]);
% 
% subplot(4,1,4);
% plot(allT,stimWave);
% xlim([allT(stimTInds(1))-width allT(stimTInds(end))+width]);

%% voltage before and after stimulation

% 
% figName = sprintf('%s_%dTrials_VoltageBeforeDuringStim%s', ...
%     experimentDescription, nTrials);
% 
% avgBefore = cell(nG,1);
% avgDuring = cell(nG,1);
% for cc = 1:Ntot
%     
%     % get the avg voltage across all trials 1 second before stim
%     vBefore = mean(allV(cc,stimTInds(1)-1000/dt:stimTInds(1)-1,:),3);
%     % get avg voltage during all trials during stim
%     vDuring = mean(allV(cc,stimTInds,:),3);
%     
%     % if stimulated population
%     if sum(cc==stimulatedNeurons)>0
%         % add to cell array
%         avgArray = avgBefore{5};
%         avgArray(end+1) = mean(vBefore);
%         avgBefore{5} = avgArray;
%         
%         avgArray = avgDuring{5};
%         avgArray(end+1) = mean(vDuring);
%         avgDuring{5} = avgArray;
%     else
%         % add to cell array of specified population
%         avgArray = avgBefore{whichpop(cc)};
%         avgArray(end+1) = mean(vBefore);
%         avgBefore{whichpop(cc)} = avgArray;
%         
%         avgArray = avgDuring{whichpop(cc)};
%         avgArray(end+1) = mean(vDuring);
%         avgDuring{whichpop(cc)} = avgArray;
%     end
% end
% 
% meanBefore = cellfun(@mean,avgBefore);
% meanDuring = cellfun(@mean,avgDuring);
% 
% % make new cell array with distributions side-by-side
% vSideBySide = cell(nG*2,1);
% inds = 1:2:nG*2;
% % make a new colormap to match this array
% cMap2 = nan(nG*2,3);
% for i=1:nG
%     vSideBySide{inds(i)} = avgBefore{i};
%     vSideBySide{inds(i)+1} = avgDuring{i};
%     cMap2(inds(i),:) = cMap(i,:);
%     cMap2(inds(i)+1,:) = cMap(i,:);
% end
% 
% meanSideBySide = cellfun(@mean,vSideBySide);
% 
% % figure out the maximum firing rate to set the y limit
% allFR = [avgBefore{:} avgDuring{:}];
% maxFR = ceil(max(allFR));
% 
% figure('Units', 'Pix', 'outerposition', [0 50 1600 1150], 'color', [1 1 1]);
% plotSpread(vSideBySide, 'distributionColors', cMap2);
% hold on;
% for i=1:nG*2
%     line([i-.4 i+.4], [meanSideBySide(i) meanSideBySide(i)], 'Color', ...
%         cMap2(i,:), 'lineWidth', 1.5);
% end
% hold off;
% newLegendCell = cell(nG*2,1);
% for i=1:nG
%     newLegendCell{inds(i)} = legendCell{i};
%     newLegendCell{inds(i)+1} = '';
% end
% legH = legend(newLegendCell);
% legH.FontSize = 16;
% a = gca;
% % a.YLim = [0 maxFR];
% a.FontSize = 14;
% ylabel('Vm (mV)', 'fontsize', 18);
% a.XTick = [];
% 
% title(strrep(figName,'_', ' '), 'Fontsize', 18,'fontweight', 'bold');
% export_fig([saveLoc, figName], '-png');




