function [Ntrials, T, dt, Ncells, Ntot, Npop, rates, times, ...
    pvtuned, p0, p2, J, r0, r2, gSyn, Esyn, taurise, taudecay, ...
    tauD, UD, tauF, UF, Fmax, Cm, gL, tau, EL, deltaT,vTpop, ...
    sigvT, vth, vre, tauref, tauw_adapt,a_adapt,b_adapt] = genParams()

Ntrials = 1;
T = 20000;
dt = 0.1;
Ncells = [400 50 25 25];
Ntot = sum(Ncells);
Npop = length(Ncells);

% curstim = 0; % no stimulation current!!

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
r0 = [2.4 0.4 0.8 0.3];
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

end