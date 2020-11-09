function Results = run_EllVerMaps(Inp, words)
%Words define the switch/case for the ellipses. Ellipse means every cell
%has the same ellipse, squished same up to being squished, rotated = same
%up to rotation, rotated and squished means both.

%Use elliptical kernels -- need global Cnns if I want hetero kernels. 

nV1 = Inp.nV1;
nThal = Inp.nThal;
nDimV1 = Inp.nDimV1;
nDimLGN = Inp.nDimLGN;
lambda1 = Inp.lambda1;
lambda2 = Inp.lambda2;
Aon = Inp.Aon;
Aoff = Inp.Aoff;

tau1 = Inp.tau(1);
tau2 = Inp.tau(2);

Cnn = globalDist(Inp.Cnn, Inp.nDimLGN,Inp.d);
Cnn_fast = globalDist(Inp.Cnn_fast, Inp.nDimLGN,Inp.d);

Gs = elongK(Inp.Xdist, Inp.sigE, Inp.sigI, words);
Gf = Inp.Gf;

sGs = svds(reshape(Gs(1,:), nDimV1, nDimV1), 1);
Gs = 1/sGs * Gs;


if length(Gf) < length(Gs)
    Gf = globalDist(Gf, Inp.nDimV1, Inp.d);
end

if length(Gf) < Inp.nV1
    sGf = svds(Gf, 1);
else
    sGf = svds(reshape(Gf(1,:),nDimV1,nDimV1),1);
end


if sGf ~= 1
    Gf = 1/sGf * Gf;
end

G1 = 0.5*Gs + tau1/(tau1+tau2)* Gf;
% G1 = Gs; 
% G2 = zeros(size(Gs)); 
G2 = tau2/(tau2+tau1)* Gs + 0.5* Gf;
K = G1 + G2;

Results.ellGs = Gs;


Cff = Cnn;
Cfn = -0.5 * Cnn;
Cnf = Cfn;

Cff_fast = Cnn_fast;
Cfn_fast = -0.5*Cnn_fast;
Cnf_fast = Cfn_fast;



%make the pure (unnormalized/unconstrained) dWdt 

% DON'T USE FFTPATH with ellipses! They differ between cortical
% cells so it doesn't work!

dWdton = @(Won, Woff) 1/4 * Aon.*((lambda1* G1) * (Won * Cnn + Woff * Cnf)+ lambda2*G2*(Won * Cnn_fast + Woff * Cnf_fast));
dWdtoff = @(Woff, Won) 1/4 * Aoff.*((lambda1* G1) * (Woff * Cff + Won * Cfn)+ lambda2*G2*(Woff * Cff_fast + Won * Cfn_fast));

%                 dWdt = @(W) 1/4 * A.*( lambda1*G1*W*p1 + lambda2*G2*W*p2 + eps*(gs*Gs + gf*Gf)*W );
lambda1fft_on = max(max(real(lambda1*fft2(G1))));
lambda2fft_on = max(max(real(lambda2*fft2(G2))));

Inp.lambda1fft_on = lambda1fft_on
Inp.lambda2fft_on = lambda2fft_on

lambda1fft_off = max(max(real(lambda1*fft2(G1))));
lambda2fft_off = max(max(real(lambda2*fft2(G2))));

Inp.lambda1fft_off = lambda1fft_off;
Inp.lambda2fft_off = lambda2fft_off;

%--------------------------------------------------------
%run Hebbian learning

if Inp.NORM
    %     while sum(sum(or(W >= Wmax, W <= 0),2)) ~= nV1*nThal % could do it
    %     this way to make all synapses inactive or
    switch Inp.IntMethod
        case 'Miller'
            [Won, Woff, WRFon, WRFoff, ts_wp, nsteps, timeUp,dt, wp] = MillerInt(Inp.Won, Inp.Woff, dWdton, dWdtoff, Inp);
            
            Inp.Nsteps = nsteps;
            Inp.dt = dt;
            
        case 'Euler'
%             [Won, wp, wpRF, ts_wp, timeUp] =  EulerMethod(Won, dWdt, Inp, dt,ts);
            [Won,Woff, WRFon,WRFoff, ts_wp, timeUp] =  EulerMethod(Inp.Won,Inp.Woff, dWdton, dWdtoff, Inp, dt, T);
%             wp1 = wp(:,:,1);
%             wp2 = wp(:,:,2);
%             wp3 = wp(:,:,3);
%             wp1RF = wpRF(:,:,1);
%             wp2RF = wpRF(:,:,2);
%             wp3RF = wpRF(:,:,3);
            Inp.dt = dt;
            
    end
Inp.endT = timeUp;

elseif ~NORM
    for t = 1:nsteps
        Won = Won +  (dt/tauW)*dWdt(Won);        
        
        if mod(t,10) == 0
            %      wtracker(:,:,t+1) = W;
            WRF = Won';
            WRF = reshape(WRF(RFinds),[Nx_arbor^d,nV1])';
            wp1(:,(t/10)+1) = WRF*e1;
%             wp2(:,(t/10)+1) = WRF*e2;
%             wp3(:,(t/10)+1) = WRF*e3;
            ts_wp((t/10)+1) = t*dt;
        end
    end
Inp.endT = t;    
end

[PrefTheta, PrefK, PrefPhi, zThetaSel, zThetaPref, zPhiSel, zPhiPref] = RFmeasure(Inp, WRFon, WRFoff);

spatialCorrRF = RFcorr(Inp, WRFon, WRFoff);

[DThetaNeighbors, DThetaRandos, DPhiNeighbors, DPhiRandos] = diffPrefTheta(Inp, zThetaPref, zPhiPref);

%Define Results Structure....................
Results.PrefTheta = PrefTheta;
Results.PrefK = PrefK;
Results.PrefPhi = PrefPhi;
Results.zThetaPref = zThetaPref;
Results.zThetaSel = zThetaSel;
Results.zPhiPref = zPhiPref;
Results.zPhiSel = zPhiSel;
Results.spatialCorrRF = spatialCorrRF;
Results.DThetaNeighbors = DThetaNeighbors;
Results.DThetaRandos = DThetaRandos;
Results.DPhiNeighbors = DPhiNeighbors;
Results.DPhiRandos = DPhiRandos;
Results.wp = wp;
Results.Won = Won;
Results.Woff = Woff;
Results.WRFon = WRFon;
Results.WRFoff = WRFoff;


[OSI, bipolarity, saturation] = OSIandBipolarity(Inp, Won, Woff, Results);
Results.OSI = OSI;
Results.bipolarity = bipolarity;
Results.saturation = saturation; 


[wp, CL] = linear2Dselectivity(Inp, Won, Woff);
[wpRF, CLRF] = ArborCovProj(Inp, Won, Woff);

Results.wp = wp;
Results.CL = CL;
Results.wpRF = wpRF;
Results.CLRF = CLRF;

end


