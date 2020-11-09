function [Won,Woff, WRFon,WRFoff, ts_wp, nsteps, timeUp,dt, wp] = MillerInt(Won,Woff, dWdton,dWdtoff, Inp, varargin)

% to create a function that will hopefully recreate the Miller 1994 paper's
% integral format. It uses a 3-step method (Birkhoff and Rota) for
% computing the integral that allows for larger time steps (lambda in the paper). 

%Formalism for computing the integral
% d/dt S(t) = lambda*F(S) --- where F(S(t)) is the final derivative after
% step 2, i.e. after normalizing somewhat. 
% S(t+1) = S(t) + lambda(23*F(t) - 16*F(t-1) + 5*F(t-2))/12 -- only works
% after 2 time steps. 

% for the first two time steps
% S(t = 1) = S(t = 0) + lambda*F(0)
% S(t= 2) = S(t=1) + lambda*(2*F(1) - F(0))

% after computation of S(t = 4) the step size and lambda were doubled, so
% that only even-numbered time-steps were computed. 
% S(t+2) = S(t) + 2*lambda(23*F(t) - 16*F(t-2) + 5*F(t-4))/12 

% Also I'm really tired of lambdas. They take too long to type. And this
% is really 1/tauW of a sort, so it should be 1/timeScale which makes me
% think of a frequency. Since it should be fancy and greek, I'll call it
% phi. 

%create time variables
tt = 4:2:200000; % Miller doubles the time step after step 4 
tt = [1,2,3,tt]; % to create the first couple of time steps
nsteps = length(tt);
ts_wp = zeros(1, nsteps);

%Pull variables from Inp structure
nV1 = Inp.nV1;
nThal = Inp.nThal;
nDimV1 = Inp.nDimV1;
nDimLGN = Inp.nDimLGN;

% p1 = Inp.p1;
% p2 = Inp.p2;
Wmax = Inp.Wmax;
% tauW = Inp.tau(3);
e1 = Inp.e1;
% e2 = Inp.e2;
% e3 = Inp.e3;
% e1RF = Inp.e1RF;
% e2RF = Inp.e2RF;
% e3RF = Inp.e3RF;
%Arbor stuff
A = Inp.A;

Aon = Inp.Aon;
Aoff = Inp.Aoff;


d = Inp.d;
RFinds = Inp.RFinds;
Nx_arbor = Inp.Nx_arbor;
R_arbor = Inp.R_arbor;
nRF = Inp.nRF;

WRFon = Won';
WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';
%make wp to be filled (nV1 by time by 3 - for slow,fast, control)
wp = zeros(nV1, length(tt), 3);
e1 = Inp.e1;
e2 = Inp.e2;
e3 = Inp.e3;
wp(:, 1, 1) = (Won -Woff)*e1;
wp(:, 1, 2) = (Won -Woff)*e2;
wp(:, 1, 3) = (Won -Woff)*e3;

%initialize variable for for loop
dt = 1; %time step for integration procedure, redefine it after first time-step. (lambda from paper)

F1on = zeros(size(Won));
F2on = zeros(size(Won));

F1off = zeros(size(Woff));
F2off = zeros(size(Woff));

if R_arbor == 5
    dt0 =  0.02;
else
    dt0 =  0.01;
end

for t = 1:nsteps
    activeOn = and(Won < Wmax.*Aon, Won.*Aon >0);%active synapses are 1, frozen = 0, syanpses freeze at Wmax*A(x-B) 
    activeOff = and(Woff < Wmax.*Aoff, Woff.*Aoff > 0);
    active = or(activeOn, activeOff); %Miller makes a distinction between
    %     synapses that are on or off being frozen.
%     activeOn = active;
%     activeOff = active; 
    
    %Step 1 Miller 1994 ----
    dW0on = dWdton(Won, Woff); % unconstrained derivative term. -- equivalent to eq 1 in Miller 94 (w/o the lambda)
%     dW0on = A.*dW0on; %don't multiply by the arbor again, already  taken into account in dWdton/off
       
    dW0on(~activeOn) = 0; 
    F1on(~activeOn) = 0;
    F2on(~activeOn) = 0;
    
    dW0off = dWdtoff(Woff, Won); % unconstrained derivative term. -- equivalent to eq 1 in Miller 94 (w/o the lambda)
%     dW0off = A.*dW0off;
       
    dW0off(~activeOff) = 0; 
    F1off(~activeOff) = 0;
    F2off(~activeOff) = 0;
    
    if t == 1
        dt = 0.01/std([dW0on(find(Aon~=0)); dW0off(find(Aoff~=0))]);
        
        if dt > dt0 %0.02
            dt = max(dt/2, dt0);
        end
        Inp.dt = dt;
        dt
    end
    
    %Step 2 Miller 1994 ----
    nActOn = sum(Aon.*activeOn,2);%for pill-box arbor this is number of active synapses onto each V1 cell
    nActOff = sum(Aoff.*activeOff,2);
%     nAct = 2*sum(A.*active,2);
    %in previous line A.*active is A_{act} of Miller-1994
    
    epsilon = sum(dW0on+dW0off,2)./(nActOn + nActOff); %epsilon will constrain the growth of W.
    epsilon = repmat(epsilon, 1, nThal);%epsilon only depends on V1-cell identity, not on LGN cell
    epsilonOn = epsilon;
    epsilonOn(~activeOn) = 0;%set epsilon-term to 0 for inactive synapses (i.e. those already on (or beyond) bounds)
    epsilonOff = epsilon;
    epsilonOff(~activeOff) = 0;
    dWon = dW0on - Aon.*epsilonOn;%from eq 2- Miller 1994 paper. It subtracts the arbor function times the constraining epsilon.
    dWoff = dW0off - Aoff.*epsilonOff; 
    
%     if t == 1
%         dt = 0.01/std([dWon(find(A)); dWoff(find(A))]);
%         
%         if dt > dt0 %0.02
%             dt = max(dt/2, dt0);
%         end
%         Inp.dt = dt;
%         dt
%     end
    
    %Step 3 Miller 1994 ----
    % what he calls F (for Final derivative - maybe) is what I've called dW. It's
    % the final derivative after step 2 at the current time step. 
    %this is the three step method derivative thing
    if t == 1
        Fon = dWon;%*dt; % current time step derivative
        W1on = Won + dt*Fon;%note that now I believe phi comes in, which is different from the words in Miller 94, but consistent with the units
        F1on = Fon; % it now becomes the previous time step derivative
        
        Foff = dWoff;
        W1off = Woff +dt*Foff;
        F1off = Foff;
    elseif t == 2
        Fon = dWon;%*dt; % current time step derivative
        W1on = Won + dt*(2*Fon - F1on);
        F2on = F1on; % two time steps previous derivative
        F1on = Fon; % previous time step derivative
        
        Foff = dWoff;%*dt; % current time step derivative
        W1off = Woff + dt*(2*Foff - F1off);
        F2off = F1off; % two time steps previous derivative
        F1off = Foff; % previous time step derivative
    else
        Fon = dWon;%*dt;
        W1on = Won + dt*(23*Fon - 16*F1on + 5*F2on)/12;
        F2on = F1on;
        F1on = Fon;
        
        Foff = dWoff;%*dt;
        W1off = Woff+ dt*(23*Foff - 16*F1off + 5*F2off)/12;
        F2off = F1off;
        F1off = Foff;
        
        if t == 4
            dt = 2*dt; %this is because the time step is now doubled -- see tt
%             dt = 2*dt; %Don't quadruple the time step. 
        end
    end
%     WAon = W1on;
%     WAoff = W1off;
    
    %Step 4 Miller 1994 ----
    %return to bound any that went strictly beyond bounds:
    W1on(W1on > Wmax.*Aon) = Wmax.*Aon(W1on > Wmax.*Aon);
    W1on(W1on < 0) = 0;
    
    W1off(W1off > Wmax.*Aoff) = Wmax.*Aoff(W1off > Wmax.*Aoff);
    W1off(W1off < 0) = 0;
    
%     if any(abs(sum(W1on+W1off,2) - 2*sum(A,2)) >= 0.0001)
%         error('Synaptic strength has changed')
%     end
    
    
    %Step 5 Miller 1994 ----
    active1On = and(W1on < Wmax.*Aon, W1on >0);%active synapses are 1, frozen = 0
    active1Off = and(W1off < Wmax.*Aoff, W1off >0);
    active1 = or(active1On, active1Off);
%     active1On = active1;
%     active1Off = active1;
    nAct1On = sum(active1On .*Aon,2);
    nAct1Off = sum(active1Off .*Aoff,2);
    nFrozOn = sum(~active1On(:)); % number of synapses that have been locked 
    nFrozOff = sum(~active1Off(:));
    
    WactOn = sum(W1on.*active1On,2);%sum of active synapses to each V1 cell
    WfrozOn = sum(W1on.*(~active1On),2);%sum of frozen synapses to each V1 cell
    
    WactOff = sum(W1off.*active1Off,2);%sum of active synapses to each V1 cell
    WfrozOff = sum(W1off.*(~active1Off),2);
    
%     gamma = repmat( (-WfrozOn-WfrozOff + nAct1On+nAct1Off)./(WactOn+WactOff) , 1, nThal);%gamma only depends on V1-cell identity; repmat( (-WfrozOn-WfrozOff + 2*sum(A,2))./(WactOn+WactOff) , 1, nThal)
%     gamma = repmat( (-WfrozOn-WfrozOff + 2*sum(A,2))./(WactOn+WactOff) , 1, nThal);
    gamma = repmat( (-WfrozOn-WfrozOff + sum(Aon,2)+sum(Aoff,2))./(WactOn+WactOff) , 1, nThal);
    gamma(gamma < 0.8) = 0.8; %see Miller step 5, end of paragraph for these limits on gamma. 
    gamma(gamma > 1.2) = 1.2;
    
    gammaOn = gamma;
    gammaOn(~active1On) = 1;%gamma-scaling is not done to frozen synapses, which are left unchanged
    gammaOff = gamma;
    gammaOff(~active1Off) = 1; %I think this should be ~activeOff because I only need to rescale Off synapses if they crossed the boundary. 
    Won = W1on.*gammaOn;%this step may push some previously active synapses beyond bounds
    Woff = W1off.*gammaOff;
    
%     Won(Won > Wmax.*A) = Wmax.*A(Won > Wmax.*A);
%     Won(Won < 0) = 0;
%     
%     Woff(Woff > Wmax.*A) = Wmax.*A(Woff > Wmax.*A);
%     Woff(Woff < 0) = 0;
    
%     W(W > Wmax.*A) = Wmax.*A(;
%     W(W < 0) = 0;'
    
%     if any(abs(sum(Won+Woff,2) - 2*sum(A,2)) >= 0.0001)
%         error('Synaptic strength has changed')
%     end
    

%     if any(any(Won > Wmax.*A))
%         pause
%     elseif any(any(Woff > Wmax.*A))
%         pause
%     end
    
    WRFon = Won';
    WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';
    WRFoff = Woff';
    WRFoff = reshape(WRFoff(RFinds),[Nx_arbor^d,nV1])';
%     wp1RF(:,t+1) = WRF*e1RF;
%     wp2RF(:,t+1) = WRF*e2RF;
%     wp3RF(:,t+1) = WRF*e3RF;
    ts_wp(t+1) = tt(t)*dt;
    
    if mod(t, 8) == 0
%         DW = sum(abs(sum(Won + Woff, 2) - 2*sum(A,2)))
        DW = sum(abs(sum(Won + Woff, 2) - sum(Aon,2) - sum(Aoff,2)))
        WRFploton = reshape(WRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
        WRFploton = permute(WRFploton, [3,1,4,2]);
        WRFploton = reshape(WRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);
        
        WRFplotoff = reshape(WRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
        WRFplotoff = permute(WRFplotoff, [3,1,4,2]);
        WRFplotoff = reshape(WRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);
    end
    
    W= Won -Woff;
    
%     wp(:,t+1,:) = W_arbored(Won, Woff, Inp);
    wp(:,t+1, 1) = W*e1;
    wp(:,t+1, 2) = W*e2;
    wp(:,t+1, 3) = W*e3;
    
    
    
%     nWon = sum(active1On(:)) 
%     nWoff = sum(active1Off(:))

%     wp2(:,t+1) = W*e2;
%     wp3(:,t+1) = W*e3;



    if sum(nAct1On+nAct1Off) <= Inp.endCondition %Stopped the simulation if 90% of neurons reach the limits 
        timeUp = t %Note this should be multiplied by dt for a true time
        if all(ts_wp(t+2:end) == 0)
            ts_wp(t+2:end) = [];
            wp(:, t+2:end, :) = [];
        else
            error('timing error')
        end
        break
    end
    
end
timeUp = t;
% WRFon = W';
% WRFon = reshape(WRFon(RFinds),[Nx_arbor^d,nV1])';


%     if t == 1
%         sdW = std(dW0(find(A))*dt)/0.01;
% %         sdW = std(dW0(find(A))*dt)/0.01
%         %tauW = sdW * tauW;
% %        tauW = min(2*sdW,1) * tauW;
%         phi = 1/sdW
%         if phi > phi0 %0.02
%             phi = max(phi/2, phi0);
%         end
%     end