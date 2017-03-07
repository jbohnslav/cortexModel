%% monosynapticExample
%  Generates four figures:
%  (E: excitatory, PV: parvalbumin, SOM: Somatostatin, VIP: vaso-intestinal
%    peptide)
%    1. E   -> E, PV, SOM, VIP
%    2. PV  -> E, PV, VIP
%    3. SOM -> E, PV, VIP
%    4. VIP -> SOM
%  Missing connections do not exist, according to this model / literature
%  This script stimulates the presynaptic neuron with current injection
%  such that the neuron spikes 10 times, waits, then spikes once more
%  (according to Litwin-Kumar et al. 2016, figure 5a)

% generate parameters
[Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams();

% overwrite some of the default model parameters
T = 2*1000; % 2 seconds long
r0 = [0 0 0 0]; % no external stimulation

% FOR FIRST FIGURE
Ncells = [2 1 1 1];
Ntot = sum(Ncells);
Npop = length(Ncells);

% generate weights
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

% CHANGE WEIGHTS
% 1 = self, 2:5 = E, PV, SOM, VIP
wstr = [0 1 1 1 1];

% set up stimulation

% other parameters
NT = round(T/dt);

stimPop = 1; % excitatory neurons
curstim = 250; % picoAmps
stimLength = 500; % ms
stimStart = 50; % ms

% stim 2
stim2Length = 50; % ms;
stim2Start = 1*1000; % ms

% big matrix, only stimulate first neuron
stimI = zeros(Ntot,NT);

stimI(1,stimStart/dt:stimStart/dt+stimLength/dt) = curstim;
stimI(1,stim2Start/dt:stim2Start/dt+stim2Length/dt) = curstim;
% JIM ADDED
allT = (dt:dt:T);
allV = nan(Ntot,NT);
allF = nan(Ntot,NT);
allD = nan(Ntot,NT);
allRise = nan(Ntot,NT);
allDecay = nan(Ntot,NT);
allExt = zeros(Ntot,NT);

% initialize other parameters

pinds = [1 1+cumsum(Ncells)]; % start index of each population (ends with Ntot+1)

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
    
    rext(pinds(pp):pinds(pp+1)-1) = r0(pp);
    
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


wext = ones(Ntot,1);
% SIMULATION
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
            
            
            % propagate spike, but only from pre-> post
            if cc==1
                % wind contains ids of downstream targets
                for kk=2:5
                    ppost = whichpop(kk);
                    xrise(pc,kk) = xrise(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                    xdecay(pc,kk) = xdecay(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                end
                for qq = 1:Npop
                    F(qq,cc) = F(qq,cc)+ UF(pc,qq)*(Fmax(pc,qq)-F(qq,cc));
                    D(qq,cc) = D(qq,cc)*UD(pc,qq);
                end
            end
        end
    end
    
    % do external excitation. NOTE as above!! This give kHz external
    % stimulation, doesn't vary with dT!!!
    %     for cc=1:Ntot
    %         while(t>nextext(cc))
    %            nextext(cc) = nextext(cc)+ exprnd(1)/rext(cc);
    %            % does this mean external excitation only goes to E neurons?
    %            xrise(1,cc) = xrise(1,cc)+wext(cc); % increment excitation
    %            xdecay(1,cc) = xdecay(1,cc)+wext(cc);
    %            allExt(cc,tt) = allExt(cc,tt)+1;
    %         end
    %     end
    
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

rates = 1000.*counts/T;

for pp=1:Npop
    fprintf('Population %d firing rate: %.2f\n', pp,...
        mean(rates(pinds(pp):pinds(pp+1)-1)));
end

times = times(1:ns);
tinds = tinds(1:ns);

% set up figure
groupNames = {'Excitatory', 'PV', 'SOM', 'VIP'};

cMap = [
    0 0.3 0.6; % blue-ish
    0 .5 0; % green
    1 .5 0; % orange
    115/255 44/255 123/255];

postYLim = [-60.2 -57.5];

saveLoc = 'D:\Analysis\cortexModel\monosynapticExample\';
figName = 'ePre';

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'Color', ...
    [1 1 1]);
% presynaptic
subplot(20,2,15:2:23);
plotH = plot(allT,allV(1,:));
a = gca;
a.XTick = [];
plotH.Color = cMap(stimPop,:);
title('Presynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{stimPop}], 'fontweight', 'bold');

% stimulation current
subplot(20,2,[25,27]);
plotH = plot(allT,stimI(1,:));
a = gca;
a.TickDir = 'out';
a.XTick = [0:500:2000];
a.YLim = [0 curstim+100];
xlabel('Time (ms)');
ylabel('Current (pA)');

% excitatory postsynaptic
subplot(20,2,[2:2:10]);
plotH = plot(allT,allV(2,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(1,:);
title('Postsynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{1}], 'fontweight', 'bold');

% PV postsynaptic
subplot(20,2,[12:2:20]);
plotH = plot(allT,allV(3,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(2,:);
ylabel([groupNames{2}], 'fontweight', 'bold');

% SOM postsynpatic
subplot(20,2,[22:2:30]);
plotH = plot(allT,allV(4,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(3,:);
ylabel([groupNames{3}], 'fontweight', 'bold');

% VIP postsynaptic
subplot(20,2,32:2:40);
plotH = plot(allT,allV(5,:));
a = gca;
a.XTick = [0:500:2000];
a.YLim = postYLim;
plotH.Color = cMap(4,:);
xlabel('Time (ms)');
ylabel([groupNames{4}], 'fontweight', 'bold');

export_fig([saveLoc, figName], '-png');

%% NOW FOR PV PRESYNAPTIC
[Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams();

% overwrite some of the default model parameters
T = 2*1000; % 2 seconds long
r0 = [0 0 0 0]; % no external stimulation

% FOR FIRST FIGURE
Ncells = [1 2 1 1];
Ntot = sum(Ncells);
Npop = length(Ncells);

% generate weights
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

% CHANGE WEIGHTS
% 1 = self, 2:5 = E, PV, SOM, VIP
wstr = [0 1 1 0 1];

% change population IDs!
whichpop = [2 1 2 3 4];
% set up stimulation

% other parameters
NT = round(T/dt);

stimPop = 2; % PV neurons
curstim = 170; % picoAmps
stimLength = 500; % ms
stimStart = 50; % ms

% stim 2
stim2Length = 50; % ms;
stim2Start = 1*1000; % ms


% big matrix, only stimulate first neuron
stimI = zeros(Ntot,NT);

stimI(1,stimStart/dt:stimStart/dt+stimLength/dt) = curstim;
stimI(1,stim2Start/dt:stim2Start/dt+stim2Length/dt) = curstim;
% JIM ADDED
allT = (dt:dt:T);
allV = nan(Ntot,NT);
allF = nan(Ntot,NT);
allD = nan(Ntot,NT);
allRise = nan(Ntot,NT);
allDecay = nan(Ntot,NT);
allExt = zeros(Ntot,NT);

% initialize other parameters

pinds = [1 1+cumsum(Ncells)]; % start index of each population (ends with Ntot+1)

% state vectors
v = -60.*ones(Ntot,1);
vT = zeros(Ntot,1);
lastSpike = -100.*ones(Ntot,1);
% whichpop = zeros(Ntot,1);

for cc = 1:Ntot
    pc = whichpop(cc);
    vT(cc) = vTpop(pp) + sigvT(pp)*randn(1,1);
    rext(cc) = r0(pp);
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


wext = ones(Ntot,1);
% SIMULATION
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
            
            
            % propagate spike, but only from pre-> post
            if cc==1
                % wind contains ids of downstream targets
                for kk=2:5
                    ppost = whichpop(kk);
                    xrise(pc,kk) = xrise(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                    xdecay(pc,kk) = xdecay(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                end
                for qq = 1:Npop
                    F(qq,cc) = F(qq,cc)+ UF(pc,qq)*(Fmax(pc,qq)-F(qq,cc));
                    D(qq,cc) = D(qq,cc)*UD(pc,qq);
                end
            end
        end
    end
    
    % do external excitation. NOTE as above!! This give kHz external
    % stimulation, doesn't vary with dT!!!
    %     for cc=1:Ntot
    %         while(t>nextext(cc))
    %            nextext(cc) = nextext(cc)+ exprnd(1)/rext(cc);
    %            % does this mean external excitation only goes to E neurons?
    %            xrise(1,cc) = xrise(1,cc)+wext(cc); % increment excitation
    %            xdecay(1,cc) = xdecay(1,cc)+wext(cc);
    %            allExt(cc,tt) = allExt(cc,tt)+1;
    %         end
    %     end
    
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

rates = 1000.*counts/T;

for pp=1:Npop
    fprintf('Population %d firing rate: %.2f\n', pp,...
        mean(rates(pinds(pp):pinds(pp+1)-1)));
end

times = times(1:ns);
tinds = tinds(1:ns);

% set up figure
groupNames = {'Excitatory', 'PV', 'SOM', 'VIP'};

cMap = [
    0 0.3 0.6; % blue-ish
    0 .5 0; % green
    1 .5 0; % orange
    115/255 44/255 123/255];

postYLim = [-65 -59];

figName = 'pvPre';

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'Color', ...
    [1 1 1]);
% presynaptic
subplot(20,2,15:2:23);
plotH = plot(allT,allV(1,:));
a = gca;
a.XTick = [];
a.YLim = [-80 20];
plotH.Color = cMap(stimPop,:);
title('Presynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{stimPop}], 'fontweight', 'bold');

% stimulation current
subplot(20,2,[25,27]);
plotH = plot(allT,stimI(1,:));
a = gca;
a.TickDir = 'out';
a.XTick = [0:500:2000];
a.YLim = [0 curstim+100];
xlabel('Time (ms)');
ylabel('Current (pA)');

% excitatory postsynaptic
subplot(20,2,[2:2:10]);
plotH = plot(allT,allV(2,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(1,:);
title('Postsynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{1}], 'fontweight', 'bold');

% PV postsynaptic
subplot(20,2,[12:2:20]);
plotH = plot(allT,allV(3,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(2,:);
ylabel([groupNames{2}], 'fontweight', 'bold');

% SOM postsynpatic
subplot(20,2,[22:2:30]);
plotH = plot(allT,allV(4,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(3,:);
ylabel([groupNames{3}], 'fontweight', 'bold');

% VIP postsynaptic
subplot(20,2,32:2:40);
plotH = plot(allT,allV(5,:));
a = gca;
a.XTick = [0:500:2000];
a.YLim = postYLim;
plotH.Color = cMap(4,:);
xlabel('Time (ms)');
ylabel([groupNames{4}], 'fontweight', 'bold');
export_fig([saveLoc, figName], '-png');
%
%% NOW FOR SOM PRESYNAPTIC
[Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams();

% overwrite some of the default model parameters
T = 2*1000; % 2 seconds long
r0 = [0 0 0 0]; % no external stimulation

% FOR FIRST FIGURE
Ncells = [1 1 2 1];
Ntot = sum(Ncells);
Npop = length(Ncells);

% generate weights
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

% CHANGE WEIGHTS
% 1 = self, 2:5 = E, PV, SOM, VIP
wstr = [0 1 1 0 1];

% change population IDs!
whichpop = [3 1 2 3 4];
% set up stimulation

% other parameters
NT = round(T/dt);

stimPop = 3; % SOM neurons
curstim = 150; % picoAmps
stimLength = 500; % ms
stimStart = 50; % ms

% stim 2
stim2Length = 50; % ms;
stim2Start = 1*1000; % ms


% big matrix, only stimulate first neuron
stimI = zeros(Ntot,NT);

stimI(1,stimStart/dt:stimStart/dt+stimLength/dt) = curstim;
stimI(1,stim2Start/dt:stim2Start/dt+stim2Length/dt) = curstim;
% JIM ADDED
allT = (dt:dt:T);
allV = nan(Ntot,NT);
allF = nan(Ntot,NT);
allD = nan(Ntot,NT);
allRise = nan(Ntot,NT);
allDecay = nan(Ntot,NT);
allExt = zeros(Ntot,NT);

% initialize other parameters

pinds = [1 1+cumsum(Ncells)]; % start index of each population (ends with Ntot+1)

% state vectors
v = -60.*ones(Ntot,1);
vT = zeros(Ntot,1);
lastSpike = -100.*ones(Ntot,1);
% whichpop = zeros(Ntot,1);

for cc = 1:Ntot
    pc = whichpop(cc);
    vT(cc) = vTpop(pp) + sigvT(pp)*randn(1,1);
    rext(cc) = r0(pp);
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


wext = ones(Ntot,1);
% SIMULATION
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
            
            
            % propagate spike, but only from pre-> post
            if cc==1
                % wind contains ids of downstream targets
                for kk=2:5
                    ppost = whichpop(kk);
                    xrise(pc,kk) = xrise(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                    xdecay(pc,kk) = xdecay(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                end
                for qq = 1:Npop
                    F(qq,cc) = F(qq,cc)+ UF(pc,qq)*(Fmax(pc,qq)-F(qq,cc));
                    D(qq,cc) = D(qq,cc)*UD(pc,qq);
                end
            end
        end
    end
    
    % do external excitation. NOTE as above!! This give kHz external
    % stimulation, doesn't vary with dT!!!
    %     for cc=1:Ntot
    %         while(t>nextext(cc))
    %            nextext(cc) = nextext(cc)+ exprnd(1)/rext(cc);
    %            % does this mean external excitation only goes to E neurons?
    %            xrise(1,cc) = xrise(1,cc)+wext(cc); % increment excitation
    %            xdecay(1,cc) = xdecay(1,cc)+wext(cc);
    %            allExt(cc,tt) = allExt(cc,tt)+1;
    %         end
    %     end
    
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

rates = 1000.*counts/T;

for pp=1:Npop
    fprintf('Population %d firing rate: %.2f\n', pp,...
        mean(rates(pinds(pp):pinds(pp+1)-1)));
end

times = times(1:ns);
tinds = tinds(1:ns);

% set up figure
groupNames = {'Excitatory', 'PV', 'SOM', 'VIP'};

cMap = [
    0 0.3 0.6; % blue-ish
    0 .5 0; % green
    1 .5 0; % orange
    115/255 44/255 123/255];

postYLim = [-65 -59];

figName = 'somPre';

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'Color', ...
    [1 1 1]);
% presynaptic
subplot(20,2,15:2:23);
plotH = plot(allT,allV(1,:));
a = gca;
a.XTick = [];
plotH.Color = cMap(stimPop,:);
title('Presynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{stimPop}], 'fontweight', 'bold');

% stimulation current
subplot(20,2,[25,27]);
plotH = plot(allT,stimI(1,:));
a = gca;
a.TickDir = 'out';
a.XTick = [0:500:2000];
a.YLim = [0 curstim+100];
xlabel('Time (ms)');
ylabel('Current (pA)');

% excitatory postsynaptic
subplot(20,2,[2:2:10]);
plotH = plot(allT,allV(2,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(1,:);
title('Postsynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{1}], 'fontweight', 'bold');

% PV postsynaptic
subplot(20,2,[12:2:20]);
plotH = plot(allT,allV(3,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(2,:);
ylabel([groupNames{2}], 'fontweight', 'bold');

% SOM postsynpatic
subplot(20,2,[22:2:30]);
plotH = plot(allT,allV(4,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(3,:);
ylabel([groupNames{3}], 'fontweight', 'bold');

% VIP postsynaptic
subplot(20,2,32:2:40);
plotH = plot(allT,allV(5,:));
a = gca;
a.XTick = [0:500:2000];
a.YLim = postYLim;
plotH.Color = cMap(4,:);
xlabel('Time (ms)');
ylabel([groupNames{4}], 'fontweight', 'bold');

export_fig([saveLoc, figName], '-png');

%% NOW FOR VIP PRESYNAPTIC
[Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams();

% overwrite some of the default model parameters
T = 2*1000; % 2 seconds long
r0 = [0 0 0 0]; % no external stimulation

% FOR FIRST FIGURE
Ncells = [1 1 1 2];
Ntot = sum(Ncells);
Npop = length(Ncells);

% generate weights
[rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2);

% CHANGE WEIGHTS
% 1 = self, 2:5 = E, PV, SOM, VIP
wstr = [0 0 0 1 0];

% change population IDs!
whichpop = [4 1 2 3 4];
% set up stimulation

% other parameters
NT = round(T/dt);

stimPop = 4; % SOM neurons
curstim = 160; % picoAmps
stimLength = 500; % ms
stimStart = 50; % ms

% stim 2
stim2Length = 50; % ms;
stim2Start = 1*1000; % ms


% big matrix, only stimulate first neuron
stimI = zeros(Ntot,NT);

stimI(1,stimStart/dt:stimStart/dt+stimLength/dt) = curstim;
stimI(1,stim2Start/dt:stim2Start/dt+stim2Length/dt) = curstim;
% JIM ADDED
allT = (dt:dt:T);
allV = nan(Ntot,NT);
allF = nan(Ntot,NT);
allD = nan(Ntot,NT);
allRise = nan(Ntot,NT);
allDecay = nan(Ntot,NT);
allExt = zeros(Ntot,NT);

% initialize other parameters

pinds = [1 1+cumsum(Ncells)]; % start index of each population (ends with Ntot+1)

% state vectors
v = -60.*ones(Ntot,1);
vT = zeros(Ntot,1);
lastSpike = -100.*ones(Ntot,1);
% whichpop = zeros(Ntot,1);

for cc = 1:Ntot
    pc = whichpop(cc);
    vT(cc) = vTpop(pp) + sigvT(pp)*randn(1,1);
    rext(cc) = r0(pp);
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


wext = ones(Ntot,1);
% SIMULATION
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
            
            
            % propagate spike, but only from pre-> post
            if cc==1
                % wind contains ids of downstream targets
                for kk=2:5
                    ppost = whichpop(kk);
                    xrise(pc,kk) = xrise(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                    xdecay(pc,kk) = xdecay(pc,kk) + ...
                        wstr(kk)*F(ppost,cc)*D(ppost,cc);
                end
                for qq = 1:Npop
                    F(qq,cc) = F(qq,cc)+ UF(pc,qq)*(Fmax(pc,qq)-F(qq,cc));
                    D(qq,cc) = D(qq,cc)*UD(pc,qq);
                end
            end
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

rates = 1000.*counts/T;

for pp=1:Npop
    fprintf('Population %d firing rate: %.2f\n', pp,...
        mean(rates(pinds(pp):pinds(pp+1)-1)));
end

times = times(1:ns);
tinds = tinds(1:ns);

% set up figure
groupNames = {'Excitatory', 'PV', 'SOM', 'VIP'};

cMap = [
    0 0.3 0.6; % blue-ish
    0 .5 0; % green
    1 .5 0; % orange
    115/255 44/255 123/255];

postYLim = [-65 -59];
figName = 'vipPre';

figure('units', 'pix', 'outerposition', [0 50 1920 1150], 'Color', ...
    [1 1 1]);
% presynaptic
subplot(20,2,15:2:23);
plotH = plot(allT,allV(1,:));
a = gca;
a.XTick = [];
plotH.Color = cMap(stimPop,:);
title('Presynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{stimPop}], 'fontweight', 'bold');

% stimulation current
subplot(20,2,[25,27]);
plotH = plot(allT,stimI(1,:));
a = gca;
a.TickDir = 'out';
a.XTick = [0:500:2000];
a.YLim = [0 curstim+100];
xlabel('Time (ms)');
ylabel('Current (pA)');

% excitatory postsynaptic
subplot(20,2,[2:2:10]);
plotH = plot(allT,allV(2,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(1,:);
title('Postsynaptic', 'fontsize', 18, 'fontweight', 'bold');
ylabel([groupNames{1}], 'fontweight', 'bold');

% PV postsynaptic
subplot(20,2,[12:2:20]);
plotH = plot(allT,allV(3,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(2,:);
ylabel([groupNames{2}], 'fontweight', 'bold');

% SOM postsynpatic
subplot(20,2,[22:2:30]);
plotH = plot(allT,allV(4,:));
a = gca;
a.XTick = [];
a.YLim = postYLim;
plotH.Color = cMap(3,:);
ylabel([groupNames{3}], 'fontweight', 'bold');

% VIP postsynaptic
subplot(20,2,32:2:40);
plotH = plot(allT,allV(5,:));
a = gca;
a.XTick = [0:500:2000];
a.YLim = postYLim;
plotH.Color = cMap(4,:);
xlabel('Time (ms)');
ylabel([groupNames{4}], 'fontweight', 'bold');
export_fig([saveLoc, figName], '-png');