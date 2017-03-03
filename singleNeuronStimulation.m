%% singleNeuronStimulation
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
%  only generate this once, average across trials using same weight matrix
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

%% set up stimulation
allT = (dt:dt:T);
% other parameters
NT = round(T/dt);
nTrials = 25;

stimPop = 4; % VIP neuron
curstim = 200; % picoAmps
stimLength = 50; % ms
stimStart = 2.5*1000; % ms

% set up groups
nStimGroups = 1;
groupToStim = 1;
groupInds = pinds(stimPop):pinds(stimPop+1)-1;

% big matrix
stimI = zeros(Ntot,NT);

idToStim = groupInds(randi(length(groupInds)));

% stimTInds = stimStart/dt:stimStart/dt+stimLength/dt;
% single square wave
% stimI(idToStim,stimTInds) = curstim;

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
figure;
plot(allT(stimTInds),squareWave);
stimI(idToStim,stimTInds) = squareWave.*curstim;
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
allSpikes = zeros(Ntot,NT);

% 100 = 10Hz imaging, 40 = 25 Hz
binSize = 40; % in milliseconds

[~,downsampledT] = downsampleSpikes(allSpikes(1,:),binSize,dt);

fRates = nan(Ntot,length(downsampledT),nTrials);

%% SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('starting sim');
for trial = 1:nTrials
    
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
        [fRate,~] = downsampleSpikes(allSpikes(cc,:),binSize,dt);
        fRates(cc,:,trial) = fRate;
    end
    
    
end

%% effect of stimulation
fprintf('Done with sim. Stimulated neuron id: %d\n', idToStim);

figure;
plot(allT./1000,mean(allV(idToStim,:,:),3));
xlabel('Time (s)');
ylabel('Vm (mV)');

%% mean response of all excitatory cells

eV = mean(allV(whichpop==1,:,:),3);

figure;
plot(allT./1000,mean(eV));

%% find large response in E on single trials
% stimTInds is the indices of stimulation

% each trial, each E neuron, get 1 second before stim
eBeforeStim = allV(whichpop==1,stimTInds(1)-1000/dt:stimTInds(1)-1,:);

eBeforeStimMean = squeeze(mean(eBeforeStim,2));
eBeforeStimStd = squeeze(std(eBeforeStim,0,2));
eDuringStim = allV(whichpop==1,stimTInds,:);
eDuringStimMean = squeeze(mean(eDuringStim,2));
eDuringStimPeak = squeeze(max(eDuringStim,[],2));

% find significant peaks on single trials by the peak voltage on that trial
% being > 2* std dev of voltage in the second before stimulation
% I think it ends up being basically a spike detector
meanDiff = eDuringStimMean-eBeforeStimMean;
peakDiff = eDuringStimPeak-eBeforeStimMean;
[maxDiff,linearInd] = max(meanDiff(:));
[maxNeuron,maxTrial] = ind2sub(size(meanDiff),linearInd);
sigDis = peakDiff>2.*eBeforeStimStd;
nSig = sum(sigDis(:));

% average pre-stim Vm across trials
eBeforeStimMeanAll = squeeze(mean(eBeforeStimMean,2));
% average stim Vm across trials
eDuringStimMeanAll = squeeze(mean(eDuringStimMean,2));
% get std dev of mean pre-stim Vm
eBeforeStimStdAll  = std(eBeforeStimMean,[],2);

sigNeurons = eDuringStimMeanAll>eDuringStimMeanAll+2.*eBeforeStimStdAll;

% get average across all trials
avgAllTrials = mean(allV(:,:,:),3);
stdAllTrials = std(allV(:,:,:),[],3);
% 
% figure;
% plot(allT./1000,allV(maxNeuron,:,maxTrial));

%%
id = 185;

figure('units', 'pix', 'outerposition', [0 50 1920 1150]);
plotH = shadedErrorBar(allT,avgAllTrials(id,:),stdAllTrials(id,:));
plotH.mainLine.LineWidth = 2;