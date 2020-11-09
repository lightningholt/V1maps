function [Won,Woff, WRFon,WRFoff, ts_wp, timeUp] =  EulerMethod(Won,Woff, dWdton,dWdtoff, Inp, dt, T, varargin)


nV1 = Inp.nV1;
nThal = Inp.nThal;
% totalSyn = nThal*nV1;
% p1 = Inp.p1;
% p2 = Inp.p2;
nsteps = round(T/dt);
if nsteps ~=Inp.Nsteps
    error('step count not right somehow')
end


ts_wp = zeros(1, nsteps);
Wmax = Inp.Wmax;
tauW = Inp.tau(3);
% e1 = Inp.e1;
% e2 = Inp.e2;
% e3 = Inp.e3;
% 
% e1RF = Inp.e1RF;
% e2RF = Inp.e2RF;
% e3RF = Inp.e3RF;

A = Inp.A;
d = Inp.d;
RFinds = Inp.RFinds;
Nx_arbor = Inp.Nx_arbor;

nRF = Inp.nRF; %sum(A(1,:)~=0);%number of presynaptic inputs to a V1 cell, also W is initialized and normalized such that sum(W,2) stays equal to nRF
% endCondition = 0.1*nRF*nV1;
%Create wp's
% wp1 = zeros(nV1, nsteps/100+1); % I don't save every times step so wp's stay smaller in size
% wp2 = wp1;
% wp3 = wp1;
% ts_wp = zeros(1, nsteps/100+1);
% wp1RF = wp1;
% wp2RF = wp1;
% wp3RF = wp1;
% % put in their initial value
% wp1(:,1) = W*e1;
% wp2(:,1) = W*e2;
% wp3(:,1) = W*e3;
% 
% WRFon = Won';
% WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';
% 
% WRF = W';
% WRF = reshape(WRF(RFinds),[Nx_arbor^d,nV1])';
% wp1RF(:,1) = WRF*e1RF;
% wp2RF(:,1) = WRF*e2RF;
% wp3RF(:,1) = WRF*e3RF;
% 

for t =1:nsteps
    
    activeOn = and(Won < Wmax.*A, Won.*A >0);%active synapses are 1, frozen = 0, syanpses freeze at Wmax*A(x-B)
    activeOff = and(Woff < Wmax.*A, Woff.*A > 0);
    active = or(activeOn, activeOff); %Miller makes a distinction between
    %     synapses that are on or off being frozen.
    %     activeOn = active;
    %     activeOff = active;
    
    %Step 1 Miller 1994 ----
    dW0on = dWdton(Won, Woff);
    dW0on(~activeOn) = 0; 
    
    dW0off = dWdtoff(Woff, Won);
    dW0off(~activeOff) = 0;
    
    nActOn = sum(A.*activeOn,2);%for pill-box arbor this is number of active synapses onto each V1 cell
    nActOff = sum(A.*activeOff,2);
    
    epsilon = sum(dW0on+dW0off,2)./(nActOn + nActOff); %epsilon will constrain the growth of W.
    epsilon = repmat(epsilon, 1, nThal);%epsilon only depends on V1-cell identity, not on LGN cell
    epsilonOn = epsilon;
    epsilonOn(~activeOn) = 0;%set epsilon-term to 0 for inactive synapses (i.e. those already on (or beyond) bounds)
    epsilonOff = epsilon;
    epsilonOff(~activeOff) = 0;
    dWon = dW0on - A.*epsilonOn;%from eq 2- Miller 1994 paper. It subtracts the arbor function times the constraining epsilon.
    dWoff = dW0off - A.*epsilonOff;
    
     %step 3 of Miller-1994 (simplified a la Euler):
    W1on = Won + dWon*dt ;%update (active) synapses
    W1off = Woff + dWoff*dt ;%update (active) synapses
    
    W1on(W1on > Wmax.*A) = Wmax.*A(W1on > Wmax.*A);
    W1on(W1on < 0) = 0;
    
    W1off(W1off > Wmax.*A) = Wmax.*A(W1off > Wmax.*A);
    W1off(W1off < 0) = 0;
    
    %Step 5 Miller 1994 ----
    active1On = and(W1on < Wmax.*A, W1on >0);%active synapses are 1, frozen = 0
    active1Off = and(W1off < Wmax.*A, W1off >0);
    active1 = or(active1On, active1Off);
%     active1On = active1;
%     active1Off = active1;
    nAct1On = sum(active1On .*A,2);
    nAct1Off = sum(active1Off .*A,2);
    nFrozOn = sum(~active1On(:)); % number of synapses that have been locked 
    nFrozOff = sum(~active1Off(:));
    
    WactOn = sum(W1on.*active1On,2);%sum of active synapses to each V1 cell
    WfrozOn = sum(W1on.*(~active1On),2);%sum of frozen synapses to each V1 cell
    
    WactOff = sum(W1off.*active1Off,2);%sum of active synapses to each V1 cell
    WfrozOff = sum(W1off.*(~active1Off),2);
    
%     gamma = repmat( (-WfrozOn-WfrozOff + nAct1On+nAct1Off)./(WactOn+WactOff) , 1, nThal);%gamma only depends on V1-cell identity; repmat( (-WfrozOn-WfrozOff + 2*sum(A,2))./(WactOn+WactOff) , 1, nThal)
    gamma = repmat( (-WfrozOn-WfrozOff + 2*sum(A,2))./(WactOn+WactOff) , 1, nThal);
    gamma(gamma < 0.8) = 0.8; %see Miller step 5, end of paragraph for these limits on gamma. 
    gamma(gamma > 1.2) = 1.2;
    
    gammaOn = gamma;
    gammaOn(~active1On) = 1;%gamma-scaling is not done to frozen synapses, which are left unchanged
    gammaOff = gamma;
    gammaOff(~active1Off) = 1; %I think this should be ~activeOff because I only need to rescale Off synapses if they crossed the boundary. 
    Won = W1on.*gammaOn;%this step may push some previously active synapses beyond bounds
    Woff = W1off.*gammaOff;
%     
%     active = and(W < Wmax, W >0);%active synapses are 1, frozen = 0
%     
%     %step 1 of Miller-1994
%     dW0 = (1/tauW)*dWdt(W);
%     dW0(~active) = 0;%set dW0 to 0 for inactive synapses (i.e. those already on (or beyond) bounds
%     
    %step 2 of Miller-1994
%     nAct = sum(A.*active,2);%for pill-box arbor this is number of active synapses onto each V1 cell
%     %in previous line A.*active is A_{act} of Miller-1994
%     
%     epsilon = sum(dW0,2)./nAct; %epsilon will constrain the growth of W.
%     epsilon = repmat(epsilon, 1, nThal);%epsilon only depends on V1-cell identity, not on LGN cell
%     epsilon(~active) = 0;%set epsilon-term to 0 for inactive synapses (i.e. those already on (or beyond) bounds)
%     dW = dW0 - A.*epsilon;%from eq 2- Miller 1994 paper. It subtracts the arbor function times the contraining epsilon.
%     
%     %step 3 of Miller-1994 (simplified a la Euler):
%     W1 = W + dW*dt ;%update (active) synapses
%     
%     %step 4 of Miller-1994:
%     %return to bound any that went strictly beyond bounds:
%     W1(W1 > Wmax) = Wmax;
%     W1(W1 < 0) = 0;
%     
%     %step 5 of Miller-1994:
%     %setting back to bounds may have changed normalization, so compensate by rescaling:
%     active1 = and(W1 < Wmax, W1 >0);%active synapses are 1, frozen = 0
%     Wact = sum(W1.*active1,2);%sum of active synapses to each V1 cell
%     Wfroz = sum(W1.*(~active1),2);%sum of frozen synapses to each V1 cell
%     gamma = repmat( (-Wfroz + nRF)./Wact , 1, nThal);%gamma only depends on V1-cell identity
%     gamma(~active1) = 1;%gamma-scaling is not done to frozen synapses, which are left unchanged
%     W = W1.*gamma;%this step may push some previously active synapses beyond bounds
    
    
    %calculate selectivities/overlaps:
    if mod(t,100) == 0
%         %      wtracker(:,:,t+1) = W;
%         
%         WRFon = Won';
%         WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';
%         wp1RF(:,(t/100)+1) = WRF*e1RF;
%         wp2RF(:,(t/100)+1) = WRF*e2RF;
%         wp3RF(:,(t/100)+1) = WRF*e3RF;
        ts_wp((t/100)+1) = t*dt;
%         
%         wp1(:,(t/100)+1) = W*e1;
%         wp2(:,(t/100)+1) = W*e2;
%         wp3(:,(t/100)+1) = W*e3;
%         
    end
    
    if sum(nAct1On+nAct1Off) <= Inp.endCondition
%     if sum(nAct) <= Inp.endCondition

        WRFon = Won';
        WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';
        WRFoff = Woff';
        WRFoff = reshape(WRFoff(RFinds),[Nx_arbor^d,nV1])';
        
        timeUp = ceil(t/100)+1; %should be multiplied by dt for a true time
%         WRF = W';
%         WRF = reshape(WRF(RFinds),[Nx_arbor^d,nV1])';
%         wp1RF(:,ceil(t/100)+1) = WRF*e1RF;
%         wp2RF(:,ceil(t/100)+1) = WRF*e2RF;
%         wp3RF(:,ceil(t/100)+1) = WRF*e3RF;
        ts_wp(ceil(t/100)+1) = t*dt;
        
%         wp1(:,ceil(t/100)+1) = W*e1;
%         wp2(:,ceil(t/100)+1) = W*e2;
%         wp3(:,ceil(t/100)+1) = W*e3;
        
        
        %{
        if all(ts_wp(ceil(t/100)+2:end) == 0)
            wp1(:, ceil(t/100)+2:end) = [];
            wp2(:, ceil(t/100)+2:end) = [];
            wp3(:, ceil(t/100)+2:end) = [];
            wp1RF(:, ceil(t/100)+2:end) = [];
            wp2RF(:, ceil(t/100)+2:end) = [];
            wp3RF(:, ceil(t/100)+2:end) = [];
            ts_wp(ceil(t/100)+2:end) = [];
            wp(:,:,1) = wp1;
            wp(:,:,2) = wp2;
            wp(:,:,3) = wp3;
            wpRF(:,:,1) = wp1RF;
            wpRF(:,:,2) = wp2RF;
            wpRF(:,:,3) = wp3RF;
        else
            error('timing error')
        end
        %}
        break
        
    end
    
    if mod(t,1000) == 0
        t
        
    end
    
end
timeUp = t;
WRFon = Won';
WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';
WRFoff = Woff';
WRFoff = reshape(WRFoff(RFinds),[Nx_arbor^d,nV1])';

% wp(:,:,1) = wp1;
% wp(:,:,2) = wp2;
% wp(:,:,3) = wp3;
% wpRF(:,:,1) = wp1RF;
% wpRF(:,:,2) = wp2RF;
%  wpRF(:,:,3) = wp3RF;
 