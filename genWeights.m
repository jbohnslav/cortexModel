function [rext,wext,wind,wipost,wstr,syncount,pinds] = ...
    genWeights(Ntot,Ncells,Npop,p0,p2,J,r0,r2)


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
pinds = [1 1+cumsum(Ncells)]; % start index of each population
for pp=1:Npop
    theta_p = linspace(-pi/2,pi/2,Ncells(pp));
    % r2 defines which inputs are tuned, and how much
    rext(pinds(pp):pinds(pp+1)-1) = r0(pp)*(1+r2(pp)*cos(2*theta_p));
end

% synaptic strength of external inputs. set to 1 and just edit gSyn
wext = ones(Ntot,1);