function Results = run_MillerVerMaps(Inp)

%under Miller's framework have the same spatial K's (G's) for both
%features; in this case Gslow = Gs. 
G1 = Inp.Gs;
G2 = Inp.Gs;
Cnn = Inp.Cnn;
Cnn_fast = Inp.Cnn_fast;


Cff = Cnn;
Cfn = -0.5 * Cnn;
Cnf = Cfn;

Cff_fast = Cnn_fast;
Cfn_fast = -0.5*Inp.Cnn_fast;
Cnf_fast = Cfn_fast;

nV1 = Inp.nV1;
nThal = Inp.nThal;
nDimV1 = Inp.nDimV1;
nDimLGN = Inp.nDimLGN;
lambda1 = Inp.lambda1;
if isfield(Inp, 'lambda3')
    lambda2 = Inp.lambda3;
else
    lambda2 = Inp.lambda2;
end

Aon = Inp.Aon;
Aoff = Inp.Aoff;

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


