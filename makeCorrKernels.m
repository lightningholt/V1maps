% make correlations and kernels

function Inp = makeCorrKernels(Inp)

ARBOR = Inp.ARBOR;
MILL94 = Inp.MILL94
d = Inp.d;
R_arbor = Inp.R_arbor;


nV1 = Inp.nV1;
nThal = Inp.nThal;
nDimV1 = Inp.nDimV1;
nDimLGN = Inp.nDimLGN;
distV1 = Inp.distV1;
gdist = Inp.gdist;


sEc = Inp.sEc;
sEc2 = Inp.sEc2;
sIc = Inp.sIc;
sIc2 = Inp.sIc2;

af = Inp.af;
aE = Inp.aE;
aI = Inp.aI;
ac = Inp.ac;

sigE = Inp.sigE;
sigI = Inp.sigI;
sigf = Inp.sigf;
a_s = Inp.a_s;

tau1 = Inp.tau(1);
tau2 = Inp.tau(2);

if Inp.d==1
    if Inp.ARBOR
        Inp.stc = 'corr1';
        p1_pars.sigEI = [sEc,sIc];
        p1_pars.aEI = [aE,aI];
        p2_pars.sigEI = [sEc2, sIc2];%0.5 * p1_pars.sigEI;
        p2_pars.aEI = p1_pars.aEI;
        [p1, p2, e1,e2, e3] = CorrProj(distV1, nThal, Inp.stc, sEc, sIc, sEc2, sIc2, aE, aI, 1);%3,5,spacing/5, spacing,  sigE, sigI, 0.5*sigE, 0.5*sigI,
        
        
        [e1RF,e2RF,e3RF] =  RFcovPCs(dx,d,p1_pars,p2_pars, R_arbor); %RF-sized vectors
        
        Cnn = p1;
        Cff = Cnn;
        Cnf = -0.5*Cnn;
        Cfn = -0.5 * Cnn;
        
        Cnn_fast = p2;
        Cff_fast = Cnn_fast;
        Cfn_fast = -0.5 * Cnn_fast;
        Cnf_fast = Cfn_fast;
    else
        Inp.stc = 'proj';
        [p1, p2, e1,e2, e3] = CorrProj(distV1, nThal, Inp.stc, sigE, sigI, 0.5*sigE, 0.5*sigI, aE, aI, 1);%3,5,spacing/5, spacing
    end
    
elseif d==2
    
    Inp.stc = 'corr2D';
    
    
    [Cnn, Cnn_fast, e1RF, e2RF, e3RF] = linearC_pc(Inp);
    Cff = Cnn;
    Cfn = -0.5 * Cnn;
    Cnf = Cfn;
    
    Cff_fast = Cnn_fast;
    Cfn_fast = -0.5 * Cnn_fast;
    Cnf_fast = Cfn_fast;
    
    CD = Cnn - Cnf;
    
    if length(Cnn)~= nV1
        Cnn_global = globalDist(Cnn, nDimLGN, Inp.d); %makeMexHat2(gdist, sEc, aE, sIc, ac, d, MILL94);%globalDist(Cnn, nDimLGN,d);
        [e1, ~] = svds(Cnn_global, 1);
    else
        Cnn_global = Cnn;
        [e1, ~] = svds(reshape(Cnn(1,:), nDimV1, nDimV1), 1);
    end
    
     if length(Cnn_fast) ~= nV1
        Cnn_fast_global =  globalDist(Cnn_fast, nDimLGN, Inp.d); %makeMexHat2(gdist, sEc2, aE, sIc2, ac, d, MILL94); %globalDist(Cnn_fast,nDimLGN, d);
        [e2, ~] = svds(Cnn_fast_global,1);
    else
        Cnn_fast_global = Cnn_fast;
        [e2, ~] = svds(reshape(Cnn_fast(1,:), nDimV1, nDimV1), 1);
     end
     
    e1 = e1/(e1'*e1)
    
    e3 = randn(length(e1),1);
    e3 = e3 - (e2'*e3)*e2 - (e1'*e3)*e1;
    e3 = e3/(e3'*e3);
    
end

Inp.Cnn = Cnn;
Inp.Cff = Cff;
Inp.Cnn_fast = Cnn_fast;
Inp.e1 = e1;
Inp.e2 = e2;
Inp.e3 = e3;
Inp.e1RF = e1RF;
Inp.e2RF = e2RF;
Inp.e3RF = e3RF;
% Inp.endCondition = 0.1*nRF*nV1; % Come back to this in a second
Inp.endCondition = 0.1*Inp.nRF; 
%--------------------------------------------------------
%--------------------------------------------------------
%make G's 

Gs = makeMexHat2(distV1,sigE,aE,sigI,aI,d, MILL94);
% Gs = globalDist(Gs, nDimV1, d);
Gs(distV1 ~= 0) = a_s*Gs(distV1~=0);
% Gs = elongK(Xdist, sigE, sigI); %makes a ellipsoid Mexican hat connection that varies from cell to cell. 
% Gs = reshape(Gs(1,:), nDimV1, nDimV1); %start simple, same kernel for all cells. 

if length(Gs) < nV1
    sGs = svds(Gs, 1);
else
    sGs = svds(reshape(Gs(1,:),nDimV1,nDimV1),1);
end


if sGs ~= 1
    Gs = 1/sGs * Gs;
end

Gf = makeMexHat2(distV1,sigf,af,0,0,d, MILL94);

if length(Gf) < length(Gs)
    Gf = globalDist(Gf, nDimV1, d);
end

if length(Gf) < nV1
    sGf = svds(Gf, 1);
else
    sGf = svds(reshape(Gf(1,:),nDimV1,nDimV1),1);
end


if sGf ~= 1
    Gf = 1/sGf * Gf;
end


% Gf = 0; 
Inp.Gs = Gs;
Inp.Gf = Gf;

% G1 = 0.5*Gs + tau1/(tau1+tau2)* Gf;
% G2 = tau2/(tau2+tau1)* Gs + 0.5* Gf;

% G1 = eye(size(Gs));
G1 = 0.5*Gs + tau1/(tau1+tau2)* Gf;
% G1 = Gs; 
% G2 = zeros(size(Gs)); 
G2 = tau2/(tau2+tau1)* Gs + 0.5* Gf;
K = G1 + G2;
% K = G1;
% CorrPlots(nDim, A, p1, p2, G1, G2)
Inp.G1 = G1;
Inp.G2 = G2;
Inp.K = K;
%--------------------------------------------------------