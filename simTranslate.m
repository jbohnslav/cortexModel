% called translate because it's my attempt to translate from Julia to
% MATLAB
% function simTranslate(T,dt,wind,wipost,wstr,wext,rext,istim,curstim)

%% GENERATE PARAMS
Ntrials = 1;
T = 20000;
dt = 0.1;
Ncells = [400 50 25 25];
Ntot = sum(Ncells);
Npop = length(Ncells);

curstim = 0; % no stimulation current!!

avcond = zeros(Npop,Ntot);
rates = zeros(Ntot,1);

times = [];
tinds = [];

pvtuned = false;

% baseline probability of connection from pre->post
p0 = [.1 .6 .6 .6;
    .6 .6 0 .6;
    .6 .6 0 .6;
    0 0 .4 0];

% tuned component of connections from pre->post (normalized to p0)
% this syntax won't work
% p_ij(dtheta) = p0_ij*(1+p2_ij*cos(2*dtheta));
p2 = [0.8 0.1 0 0;
      0.8 0.1 0 0;
      0 0 0 0;
      0 0 0 0];
  
% synaptic strength for existing connection pre-> post. Set to 1 and just
% use gsyn
J = [1 1 1 1;
     1 1 1 1;
     1 1 1 1;
     1 1 1 1];
 
% set up external input rates
r0 = [2.4 0.4 0.8 0.3]*1;
r2 = [0.2 0.2 0 0];

if ~pvtuned
    p2(2,:) = 0;
    r2(2) = 0;
end

% synaptic parameters, Matrix for pre- and post-pop (nanoSiemens)
% THIS IS TRANSPOSED FROM PAPER!
gSyn = [.1/60 .3/60 .05/60 .1/60;
        3/22 3/22 0 .6/22;
        1.5/22 1/22 0 2.5/22;
        0 0 3/22 0]*1000;
    
% reversal potential, milliVolts
Esyn = [0 0 0 0;
    -67 -67 -67 -67;
    -67 -67 -67 -67;
    -67 -67 -67 -67];

% rise time, milliseconds
taurise = [.5 .5 .5 .5;
    .5 .5 .5 .5;
    1 1 1 1;
    1 1 1 1];

% decay time, milliseconds
taudecay = [2 2 2 2;
    3 3 3 3;
    4 4 4 4;
    4 4 4 4];

% depression time constant, milliseconds
tauD = [800 800 800 800;
    800 800 800 800;
    800 800 800 800;
    800 800 800 800];

% amount of depression
UD = [.75 .75 1 1;
    .9 .9 .9 .9;
    1 1 1 1;
    1 1 1 1];

% facilitation time constant (ms)
tauF =  [200 200 200 200;
    200 200 200 200;
    200 200 200 200;
    200 200 200 200];

% amount of facilitation
UF = [0 0 .5 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];

% maximum facilitation
Fmax = [2 2 2 2;
    2 2 2 2;
    2 2 2 2;
    2 2 2 2];

% cell parameters, one for each population
Cm = [180 80 80 80]; % capacitance, pF
gL = [1/0.16 1/.1 1/.2 1/.2]; % leak conductance, nanoSiemens

tau = Cm./gL; % membrane time constant

EL = [-60 -60 -60 -60]; % leak voltage, mV
deltaT = [1 .25 1 1]; % aEIF slope parameter, mV
vTpop = [-40 -40 -45 -45]; % aEIF threshold, mV
sigvT = [3 3 3 3]; % std dev of threshold
vth = [20 20 20 20]; % maximum voltage, mV
vre = [-60 -60 -60 -60]; % reset voltage, mV
tauref = [2 2 2 2]; % refractory period, ms
tauw_adapt = [150 150 150 150]; % adaptation timescale, ms

a_adapt = [4 0 4 4]; % adaptation slope, nanoSiemens
b_adapt = 10.*[.8 0 .8 .8]; % adaptation increment, pA
%% genweights: generate recurrent and external weights and external input
%  rates

% set up recurrent weight matrix
syncount = 1;
Maxw = round(Ntot*Ntot*.3); % maximum number of weights in weight matrix
% column of w corresponding to the start of the ith neuron's projections
wind = zeros(Ntot+1, 1); 
wipost = zeros(Maxw,1);
wstr = zeros(Maxw,1);
% iterate over populations
for pp=1:Npop
    theta_p = linspace(-pi/2, pi/2, Ncells(pp));
    % iterate over cells in this population
    for cc = 1:Ncells(pp)
        wind(cc + sum(Ncells(1:pp-1))) = syncount; 
        
        theta = theta_p(cc);
        
        % iterate over all populations again
        for qq = 1:Npop
            theta_q = linspace(-pi/2,pi/2,Ncells(qq));
            
            dtheta = theta - theta_q;
            % find the minimum absolute difference between this cell's
            % orientation tuning and the orientation tuning of cells in the
            % currently iterated-over population
            dtheta = min([abs(dtheta),pi-abs(dtheta)]);
            % probability of connection between this neuron and every one
            % in the target population
            prob = p0(pp,qq)*(1+p2(pp,qq)*cos(2*dtheta));
            % find the indices of non-zero connections
            iconns = find(rand(1,Ncells(qq)) < prob)+sum(Ncells(1:qq-1));
            wipost(syncount:syncount+length(iconns)-1)= iconns;
            % set pre->post connection weights to 1!! the actual weight is
            % set as the synaptic conductance gSyn when the simulation runs
            wstr(syncount:syncount+length(iconns)-1) = J(pp,qq);
            syncount = syncount+length(iconns);
        end
    end
end
% not sure why I have to change this... indexing in julia?
syncount = syncount-1;
% clip empty synapses
wind(Ntot+1) = syncount;
wipost = wipost(1:syncount);
wstr = wstr(1:syncount);

rext = zeros(Ntot,1);
% pinds is already defined
for pp=1:Npop
    theta_p = linspace(-pi/2,pi/2,Ncells(pp));
    % r2 defines which inputs are tuned, and how much
    rext(pinds(pp):pinds(pp+1)-1) = r0(pp)*(1+r2(pp)*cos(2*theta_p));
end

% synaptic strength of external inputs. set to 1 and just edit gSyn
wext = ones(Ntot,1);

%% initialize other parameters
Ntot = sum(Ncells);
Npop = length(Ncells);
pinds = [1 1+cumsum(Ncells)]; % start index of each population (ends with Ntot+1)

% simulation parameters
NT = round(T/dt);
recstart = 2000; % time at which to start analysis, to remove transient (?)
irecstart = round(recstart/dt);

% external stimulation
curext = zeros(Ntot,1); % constant bias current to one population
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
maxns = round(1,T*Ntot*.05);
times = zeros(1,maxns);
tinds = zeros(1,maxns);
vall = zeros(1);

% average conductance each neuron received
avcond = zeros(Npop,Ntot);
%% JIM ADDED
allT = (dt:dt:T);
allV = nan(Ntot,NT);
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
            dv = curext(cc) - w_adapt(cc)-gL(pc)*(v(cc)-EL(pc))+...
                gL(pc)*deltaT(pc)*exp((v(cc)-vT(cc))/deltaT(pc));
            
            % ADD SYNAPTIC CURRENTS HERE! Current depends on population
            for qq=1:Npop
                % why is this a minus sign??
                dv = dv - gSyn(qq,pc)*(xdecay(qq,cc)-xrise(qq,cc))/(taudecay(qq,pc)-taurise(qq,pc));
            end
            
            % change current to voltage
            v(cc) = v(cc)+dt*dv/Cm(pc);
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
        end
    end
    
    if tt>irecstart
        for cc = 1:Ntot
           pc = whichpop(cc);
           for qq=1:Npop
              avcond(qq,cc) = avcond(qq,cc) + ...
                  gSyn(qq,pc)*(xdecay(qq,cc)-xrise(qq,cc))/(taudecay(qq,pc)-taurise(qq,pc)); 
           end
        end
    end
    
    % JIM ADDED
    allV(:,tt) = v;
end % end sim

disp('100%');
disp('Calculating rates');
counts = zeros(1,Ntot);
for cc=1:Ntot
    counts(cc) = sum((tinds==cc)*(times>recstart));
end

rates = 1000.*counts/(T-recstart);

for pp=1:Npop
    fprintf('Population %d firing rate: %.2f\n', pp,...
        mean(rates(pinds(pp):pinds(pp+1)-1)));
end

times = times(1:ns);
tinds = tinds(1:ns);
avcond = avcond/(NT-irecstart);

