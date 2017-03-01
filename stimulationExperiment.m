%% stimulationExperiment
close all; clear all; clc

%% GENERATE PARAMS
[Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams();
%% genweights: generate recurrent and external weights and external input
%  rates
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

%% set up stimulation

% other parameters
NT = round(T/dt);

stimPop = 1; % excitatory neurons
curstim = 1000; % picoAmps
stimLength = 5; % ms
stimStart = 5*1000; % ms

% set up groups
nGroups = 4;
groupToStim = 1;
groupInds = pinds(stimPop):pinds(stimPop+1)-1;

% big matrix
stimI = zeros(Ntot,NT);

groupIndsRand = groupInds(randperm(length(groupInds)));

nPerGroup = length(groupInds)/nGroups;
for i=1:nGroups
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

% simulation parameters

% external stimulation
curext = zeros(Ntot,1); % constant bias current to one population
istim = 1; % index identity of population to be stimulated, in this case E
curext(pinds(istim):pinds(istim+1)-1) = curstim;

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
%% SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('starting sim');
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
    
%     if tt>irecstart
%         for cc = 1:Ntot
%            pc = whichpop(cc);
%            for qq=1:Npop
%               avcond(qq,cc) = avcond(qq,cc) + ...
%                   gSyn(qq,pc)*(xdecay(qq,cc)-xrise(qq,cc))/(taudecay(qq,pc)-taurise(qq,pc)); 
%            end
%         end
%     end
    
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
%% raster figure

figure;
subplot(4,3,1:9);
for i=1:Ntot
    spikes = allSpikes(i,:);
    spikeTimes = allT(find(spikes));
    if ~isempty(spikeTimes)
        if sum(stimulatedNeurons==i)
            for j=1:length(spikeTimes)
                line([spikeTimes(j) spikeTimes(j)], [i-1 i], 'color', ...
                    'b');
            end
        else
            for j=1:length(spikeTimes)
                line([spikeTimes(j) spikeTimes(j)], [i-1 i], 'color', ...
                    'k');
            end
        end
    end
end
ylim([0 Ntot]);
subplot(4,3,10:12);
[row,col] = find(stimI==max(max(stimI,[],2)),1);
plot(allT,stimI(row,:));

%% average firing rate figure
binSize = 100; % in milliseconds
groupNames = {'E', 'PV', 'SOM', 'VIP', 'stimE'};
[~,downsampledT] = downsampleSpikes(allSpikes(1,:),binSize,dt);

fRates = nan(Ntot,length(downsampledT));
disp('Binning firing rates');
for cc=1:Ntot
    [fRate,~] = downsampleSpikes(allSpikes(cc,:),binSize,dt);
    fRates(cc,:) = fRate;
end

groupAverages = zeros(nGroups+1,length(downsampledT));
gNums = zeros(nGroups+1,1);
gNum = 1;
for i=1:Ntot
    if sum(stimulatedNeurons==i)>0
        groupAverages(nGroups+1,:) = groupAverages(nGroups+1,:)+...
            fRates(i,:);
        gNums(nGroups+1) = gNums(nGroups+1)+1;
    else
        groupAverages(whichpop(i),:) = groupAverages(whichpop(i),:)+...
            fRates(i,:);
        gNums(whichpop(i)) = gNums(whichpop(i))+1;
    end
end
% firing rate averages, in hertz
groupAverages = groupAverages./gNums.*binSize*dt;


figure;
subplot(4,3,1:9);
hold on;
for i=1:nGroups+1
    plot(downsampledT,groupAverages(i,:));
end
legend(groupNames);
