% Analysis functions:
% RFcorr -- finds correlations between pairs out to distance R
% diffPrefTheta -- Finds differences between neigbhors preferred oris
% absolut_Phase_calc -- find absolute phase
% absoluteSpatialCorr -- finds absolute Spatial correlation
% correlationLengthCalc -- finds correlation length


%Figure functions:
% linearFeatureMaps -- makes good figure for linear selectivities
% comp_Results -- makes good general figures
% OriPref_comp -- makes good ori figures
% higherOrderPC_slecivities -- makes higher order PC selectivity maps
% RFcompMillandMyframe -- compares RFs between my framework and miller
% ExAbs_RelCorr -- makes examples of RFs which show absolute and relative
%                 correlation schematically plus histogrmas of Correlation
% RF_plots -- makes full field RF plot
% phasePref_AbsVsRel -- compares maps of absolute and relative phase
% visualizeCandK -- Makes visualizations of 2D C and K



%I'm making this to make better figures for my dissertation
%2 - good lambda2 = 10
%3 - good 
%4 - linCfast not great
%5 - lambda2 7 not great Cfast
%6 - lambda2 = 8.5 not great Cfast 
%7 - lambda2 = 9 better than 8.5

% 10 - lambda2 =10 for 64x64 neurons, lam3= 1
% 11 - lambda2 = 10 for 32x32 neurons, lam3 = 10 -- makes very smooth
%      linear map for C_fast, bad looking RFs
% 12 - lambda2 = 7 for 32x32 neurons, lam3 = 7 -- makes things very similar
%      to 11, but the linear selectivity for the fast feature in dyn looks much
%      weaker.

% 13 lambda2= 7 for both. same result hold
% 14 lambda2=10 
% 15 lambda2= 10 -- Fixed eRFs. Now both features were S&P looking,
%    although I was looking at the selectivity to e's.. not eRFs/,
% 16 lambda2 = 5, still looked S&Pery -- this is wpRF, wp looks pretty
%    good, long range correlation in smooth
% 17 lambda2 = 1 strong smooth map, weak fast map
% 18 lambda2 = 3 strong smooth map, better fast map
% 19 - lambda2 = 7 -- strengths looked pretty good, but the smooth map
%     wasn't very smooth
% 20 - lambda2 = 6 -- looks better than 7, but still not smooth

% 21 - lambda2 = 5.5, lambda3 = 2.5 -- goldilocks


function [Inp, Big_Results] = C10VisMaps(ii, lambda2, lambda3, dirname)

% clear
% ii = 7;
rc = 0.20; 
rI = 0.3; 
c_a = 0.5; %this is the factor that scales the dendritic radius

% lambda2 = 10%.24; %25/3.5;
lambda1 = 1; %lambda2 is given as first input to function. (normally 10);

Inp = gridDistances(rc, rI, c_a);

tau1 = 100;
tau2 = 1;
gs = 1/tau1; %this is gamma_s from Yashar's notes -- arbitrarily setting gamma_s = 1/tau1
gf = 1/tau2; % this is gamma_f from Yashar's notes

dt = 0.1*tau2; % dt sets the time scale -- currently 1/10 of tau2.
tauW = 1000;
scaleTauW = tauW/(lambda1*Inp.aE);
T = round(5*10*tauW/(lambda1*Inp.aE));


Inp.tau = [tau1; tau2; tauW]; %tau1, tau2, tauW defined in ConDefinitions.mat
Inp.scaleTau = tauW/(lambda1*Inp.aE);
Inp.lambda1 = lambda1; 
Inp.lambda2 = lambda2;
Inp.lambda3 = lambda3;

switch Inp.IntMethod
    case 'Euler'
        Inp.T = T;
        Inp.Nsteps = round(T/dt);
end
%make arbor function



Inp = makeCorrKernels(Inp);

nV1 = Inp.nV1;
nThal = Inp.nThal;
nDimV1= Inp.nDimV1;
nDimLGN = Inp.nDimLGN;

Aon = Inp.Aon;
Aoff= Inp.Aoff;

G1 = Inp.G1;
G2 = Inp.G2;

Cnn = Inp.Cnn;
Cff = Cnn;
Cnf = -0.5 * Cnn;
Cfn = Cnf;

Cnn_fast = Inp.Cnn_fast;
Cff_fast = Cnn_fast;
Cfn_fast = -0.5 * Cnn_fast;
Cnf_fast = Cfn_fast;


%make the pure (unnormalized/unconstrained) dWdt 
switch Inp.stc
    case {'proj','corr1'}
        switch Inp.fftpath
            case 'brute-force'
                dWdt = @(W) 1/4 * A.*( lambda1*G1*W*p1 + lambda2*G2*W*p2); %this adds a more complete term. ;
                lambda1fft = max(real(lambda1*fft(G1(1,:)')))
                lambda2fft = max(real(lambda2*fft(G2(1,:)')))
                Inp.lambda1fft = lambda1fft;
                Inp.lambda2fft = lambda2fft;
                
            case 'yes' %currently this case assumes eps is 0
                if d~=1
                    error('d must be 1, in case corr1');
                end
                GCff_fft = fft(G1(1,:)') * fft(p1(1,:)); %nV1 x nThal matrix
                GCnf_fft = fft(G2(1,:)') * fft(p2(1,:)); %nV1 x nThal matrix
                dWdt = @(W) 1/4 * fastdWdt2(W,A, lambda1*GCff_fft + lambda2*GCnf_fft);
                lambda1fft = max(real(lambda1*GCff_fft(:)))
                lambda2fft = max(real(lambda2*GCnf_fft(:)))
                Inp.lambda1fft = lambda1fft;
                Inp.lambda2fft = lambda2fft;
        end
    case 'corr2D' %currently this case assumes eps is 0
        switch Inp.fftpath
            case 'yes'
%                 on cells
%                 G1fft = reshape(G1(:,1), nDimV1, nDimV1); %gets first matrix describing connectivity between (-L, -L) neuron and rest of network in V1
                
                if length(G1) ~= nV1
                    G1fft = reshape(fft2(G1), nV1, 1); %takes the fft in 2d and then vectorizes the results
                    G2fft = reshape(fft2(G2), nV1, 1);
                    
                    Cnnfft = reshape(fft2(Cnn),nThal, 1); %takes the fft in 2d and then vectorizes the results
                    fastCnn_fft = reshape(fft2(Cnn_fast), nThal,1);
                    %                 Cnffft = reshape(Cnf(:,1), nDimLGN, nDimLGN); %gets first matrix describing connectivity between (-L, -L) neuron and the rest of network in Thalamus
                    Cnffft = reshape(fft2(Cnf),nThal, 1);
                    fastCnf_fft = reshape(fft2(Cnf_fast), nThal,1);
                    
                    GCnn_fft = kron(Cnnfft, G1fft); %takes the (dirac? tensor?) product of the two. So now dimensionality is (nV1*nThal)x1
                    fastGCnn_fft = kron(fastCnn_fft, G2fft);
                    
                    
                    GCnf_fft = kron(Cnffft, G1fft);
                    
                    fastGCnf_fft = kron(fastCnf_fft, G2fft);
                    
                    Cff_fft = reshape(fft2(Cff),nThal, 1); %takes the fft in 2d and then vectorizes the result;
                    fastCff_fft = reshape(fft2(Cff_fast), nThal,1);
                    %                 Cfn_fft = reshape(Cfn(:,1), nDimLGN, nDimLGN); %gets first matrix describing connectivity between (-L, -L) neuron and the rest of network in Thalamus
                    Cfn_fft = reshape(fft2(Cfn),nThal, 1);
                    fastCfn_fft = reshape(fft2(Cfn_fast), nThal, 1);
                    
                    GCff_fft = kron(Cff_fft, G1fft); %takes the (dirac? tensor?) product of the two. So now dimensionality is (nV1*nThal)x1
                    fastGCff_fft = kron(fastCff_fft, G2fft);
                    
                    GCfn_fft = kron(Cfn_fft, G1fft);
                    fastGCfn_fft = kron(fastCfn_fft, G2fft);
                    
                else
                    
                    G1fft = fft2(reshape(G1, [nDimV1, nDimV1, nDimV1, nDimV1]));
                    G2fft = fft2(reshape(G2, [nDimV1, nDimV1, nDimV1, nDimV1]));
                    Cnnfft = fft2(reshape(Cnn_global, [nDimLGN, nDimLGN, nDimLGN, nDimLGN]));
                    Cnffft = -0.5 * Cnnfft;
                    fastCnn_fft = fft2(reshape(Cnn_fast_global, [nDimLGN, nDimLGN, nDimLGN, nDimLGN]));
                    fastCnf_fft = -0.5 * fastCnn_fft;
                    
                    GCnn_fft = Cnnfft .* G1fft;
                    fastGCnn_fft = fastCnn_fft .* G2fft;
                    
                    GCnf_fft = Cnffft .* G1fft;
                    fastGCnf_fft = fastCnf_fft .* G2fft;
                    
                    Cff_fft = fft2(reshape(Cnn_global, [nDimLGN, nDimLGN, nDimLGN, nDimLGN]));
                    fastCff_fft = fft2(reshape(Cnn_global, [nDimLGN, nDimLGN, nDimLGN, nDimLGN]));
                    Cfn_fft = -0.5 * Cff_fft;
                    fastCfn_fft = -0.5 * fastCff_fft;
                    
                    GCff_fft = Cff_fft .* G1fft;
                    fastGCff_fft = fastCff_fft .* G2fft;
                    
                    GCfn_fft = Cfn_fft .* G1fft;
                    fastGCfn_fft = fastCfn_fft .* G2fft;
                    
%                     dWdton = @(Won, Woff) 1/4 *(fastdWdt2(Won,Aon, (1-MvA)*(lambda1*GCnn_fft+lambda2*fastGCnn_fft) + MvA*PCnn_fft) + fastdWdt2(Woff,Aon, (1-MvA)*(lambda1*GCnf_fft+lambda2*fastGCnf_fft) +MvA*PCnf_fft));
%                     dWdtoff = @(Woff, Won) 1/4 *(fastdWdt2(Woff,Aoff, (1-MvA)*(lambda1*GCff_fft+lambda2*fastGCff_fft) + MvA*PCff_fft) + fastdWdt2(Won,Aoff, (1-MvA)*(lambda1*GCfn_fft+lambda2*fastGCfn_fft) +MvA*PCfn_fft));
%                     
%                     dWdton = @(Won, Woff) ((1/4)* Aon .* ifft2(((1-MvA) * (lambda1 * GCnn_fft + lambda2 * fastGCnn_fft) + MvA*PCnn_fft) .* fft2(Won) + ((1-MvA) * (lambda1 * GCnf_fft + lambda2 * fastGCnf_fft) + MvA*PCnf_fft).* fft2(Woff)));
%                     dWdtoff = @(Woff, Won) ((1/4)* Aoff .* ifft2(((1-MvA) * (lambda1 * GCff_fft + lambda2 * fastGCff_fft) + MvA*PCff_fft) .* fft2(Woff) + ((1-MvA) * (lambda1 * GCfn_fft + lambda2 * fastGCfn_fft) + MvA*PCfn_fft).* fft2(Won)));
                end
                
                GCnn_fft = reshape(GCnn_fft,[nDimV1,nDimV1,nDimLGN,nDimLGN]); %think of this as like 2 2-vector indices. In the 1D case, those 2-vectors become 1-vectors.
                GCnf_fft = reshape(GCnf_fft, [nDimV1,nDimV1,nDimLGN,nDimLGN]);
                fastGCnn_fft = reshape(fastGCnn_fft, [nDimV1,nDimV1,nDimLGN,nDimLGN]);
                fastGCnf_fft = reshape(fastGCnf_fft, [nDimV1,nDimV1,nDimLGN,nDimLGN]);
                
                GCff_fft = reshape(GCff_fft,[nDimV1,nDimV1,nDimLGN,nDimLGN]); %think of this as like 2 2-vector indices. In the 1D case, those 2-vectors become 1-vectors.
                fastGCff_fft = reshape(fastGCff_fft,[nDimV1,nDimV1,nDimLGN,nDimLGN]);
                GCfn_fft = reshape(GCfn_fft, [nDimV1,nDimV1,nDimLGN,nDimLGN]);
                fastGCfn_fft = reshape(fastGCfn_fft,[nDimV1,nDimV1,nDimLGN,nDimLGN]);
                
                dWdtoff = @(Woff, Won) 1/4 *(fastdWdt2(Woff,Aoff, lambda1*GCff_fft+lambda2*fastGCff_fft , nDimV1, nDimLGN) + fastdWdt2(Won,Aoff, lambda1*GCfn_fft+lambda2*fastGCfn_fft, nDimV1, nDimLGN));
                dWdton = @(Won, Woff) 1/4 *(fastdWdt2(Won, Aon, lambda1*GCnn_fft+lambda2*fastGCnn_fft, nDimV1, nDimLGN) + fastdWdt2(Woff, Aon, lambda1*GCnf_fft+lambda2*fastGCnf_fft, nDimV1, nDimLGN));
                
                lambda1fft_on = max(real(lambda1*GCnn_fft(:)))
                lambda2fft_on = max(real(lambda2*fastGCnn_fft(:)))
                
                lambda1Diff = max(real(lambda1*(GCnn_fft(:) - GCnf_fft(:))))
                lambda2Diff = max(real(lambda2*(fastGCnn_fft(:) - fastGCnf_fft(:))))
                
                lambda1fft_off = max(real(lambda1*GCff_fft(:)))
                lambda2fft_off = max(real(lambda2*fastGCff_fft(:)))
                
                Inp.lambda1fft_on = lambda1fft_on;
                Inp.lambda2fft_on = lambda2fft_on;
                Inp.lambda1fft_off = lambda1fft_off;
                Inp.lambda2fft_off = lambda2fft_off;
                
            case 'brute-force'
                 
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
        end
        
        
end

% Inp.lambda1fft = lambda1fft;
% Inp.lambda2fft = lambda2fft;
%--------------------------------------------------------
%--------------------------------------------------------
%set the initial condition for W
Wmax = 4;
Wmean = 1;
Wnoise = 0.2;
% sigW = 0.2*dxConv;d %*96; %this is from Antolik and Bednar paper

Won = Wmean + Wnoise*(2*rand(nV1, nThal)-1);%Uniform distribution on [Wmean-Wnoise, Wmean+Wnoise];
% Wbase = exp(-gdist.^2./(2*sigW^2));
% Won = rand(nV1, nThal).*Wbase; %Consistent with Antolik and Bednar

Won = Aon.*Won;

Woff = Wmean + Wnoise*(2*rand(nV1, nThal)-1);%Uniform distribution on [Wmean-Wnoise, Wmean+Wnoise];
% Woff = rand(nV1, nThal).*Wbase; %consistent with Antolik and Bednar
Woff = Aoff.*Woff;

if Inp.NORM %rescale to set the sum of weights to V1 cells equal to sum(A,2)
%     W1on = repmat(2*sum(A,2)./sum(Won+Woff,2) ,1, nThal).*Won; %see Miller Algorithm step 5. 
%     Woff = repmat(2*sum(A,2)./sum(Won+Woff,2),1,nThal).*Woff; %the 2 is cause there are now two grids of LGN
    W1on = repmat(sum(Aon+Aoff,2)./sum(Won+Woff,2) ,1, nThal).*Won; %see Miller Algorithm step 5. 
    W1off = repmat(sum(Aon+Aoff,2)./sum(Won+Woff,2),1,nThal).*Woff; %the 2 is cause there are now two grids of LGN
    
    Won = W1on;
    Woff = W1off;
end

Inp.Won = Won;
Inp.Woff = Woff;
Inp.Wmax = Wmax;
%--------------------------------------------------------
%run Hebbian learning

if Inp.NORM
    %     while sum(sum(or(W >= Wmax, W <= 0),2)) ~= nV1*nThal % could do it
    %     this way to make all synapses inactive or
    switch Inp.IntMethod
        case 'Miller'
            [Won, Woff, WRFon, WRFoff, ts_wp, nsteps, timeUp,dt, wp] = MillerInt(Won, Woff, dWdton, dWdtoff, Inp);
            
            Inp.Nsteps = nsteps;
            Inp.dt = dt;
            
        case 'Euler'
%             [Won, wp, wpRF, ts_wp, timeUp] =  EulerMethod(Won, dWdt, Inp, dt,ts);
            [Won,Woff, WRFon,WRFoff, ts_wp, timeUp] =  EulerMethod(Won,Woff, dWdton, dWdtoff, Inp, dt, T);
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

%========================================================================
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

[OSI, bipolarity, saturation] = OSIandBipolarity(Inp, Won, Woff, Results);
Results.OSI = OSI;
Results.bipolarity = bipolarity;
Results.saturation = saturation; 

if Inp.d == 1
    figname = 'Feedback_1Dring';
elseif Inp.d ==2
    figname = 'MyFrame_2Dgrid';
end

% if ARBOR 
%     figname = [figname,'_Arbor'];
% end
if Inp.NORM 
    figname = [figname,'_Norm'];
end

switch Inp.fftpath 
    case 'yes'
    figname =[figname, '_fft'];
end
switch Inp.IntMethod 
    case 'Miller'
    figname = [figname, '_Miller'];
    case 'Euler'
    figname = [figname, '_Euler'];
end


figname = [figname, num2str(ii)];

[wp, CL] = linear2Dselectivity(Inp, Won, Woff);
[wpRF, CLRF] = ArborCovProj(Inp, Won, Woff);

Results.wp = wp;
Results.CL = CL;
Results.wpRF = wpRF;
Results.CLRF = CLRF;
Results.Won = Won;
Results.Woff = Woff;
Results.WRFon = WRFon;
Results.WRFoff = WRFoff;

Mill_Results = run_MillerVerMaps(Inp);
% ell_Results = run_EllVerMaps(Inp, 'ellipse');
% rot_Results = run_EllVerMaps(Inp, 'rotated');
% squi_Results = run_EllVerMaps(Inp, 'squished');
% squi_rot_Results = run_EllVerMaps(Inp, 'rotated and squished');

Mill_Results.absSpatCorr = absoluteSpatialCorr(Mill_Results, Inp);
Results.absSpatCorr = absoluteSpatialCorr(Results, Inp);

Big_Results.Res= Results;
Big_Results.Mill_Res = Mill_Results;
% Big_Results.ell_Res = ell_Results;
% Big_Results.rot_Res = rot_Results;
% Big_Results.squi_Res = squi_Results;
% Big_Results.squi_rot_Res = squi_rot_Results;

save(['data/',figname], 'Inp', 'Big_Results', '-v7.3');
comp_Results(Big_Results, Inp, ii, dirname);


makePDFs(figname,Inp, WRFon, WRFoff, Results) %makes the standard pdf figure that I like.

