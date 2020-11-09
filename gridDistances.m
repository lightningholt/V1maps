% Makes distances for the grid/ring
function Inp = gridDistances(rc, rI, c_a)

Inp.rc = rc;
Inp.rI = rI;

af = 1; %0.67;%0.3

NORM = 1; %Doesn't normalize W (Norm=0), or uses Miller norm (NORM = 1);
ARBOR = 1; %uses the arbor function or doesn't (ARBOR = 0). 
MILL94 = 1;
PBC= 1; % Periodic Boundary Conditions
IntMethod = 'Miller'; %My previous results were using the Euler method. 
fftpath = 'yes'; %Switch between fft ('yes') and not ('brute-force').  -- can't use this with heterogeneous Ks

d = 2;%topographic/retinotopy dimension 

if d==1
    L = 10
    nDim = 1024;%256; %Miller 1994 uses 1024 neurons in a 32x32 grid
    R_arbor = 1; %sets length of arbor
   
elseif d==2
    nDim = 32; % want this to be 16
    L = 16 %Antolik length = 1 unit/side; %nDim/2 %0.3125 %To make spacing the same for 1D case and 2D case
    R_arbor = 6*(2*L/nDim);% 4 * dx -- Gives me a diameter of 11, Miller has an effective diameter of 11.3. 
%     R_arbor = 0.2*(2*L/nDim)/(1/96); %19 %from Antolik Paper
end

D_arbor = 2*floor(R_arbor/(2*L/nDim))+1; %11; 
nDimV1 = nDim;
nDimLGN = nDim; %But there are two grids of LGN cells (ON and OFF center cells)
dx = 2*L/(nDim);
nV1 = nDim^d; %number of V1 neurons -- one dimensional currently
nThal = nV1; %sets size of input layer. =nV1 means equal number of input and output neurons

aE = 1; %this is one in Miller 94
aI = 1/9;

a_s = 0.5; % multiplies the G's (interaction term) by something less than 1 for distances greater than 0;
ac = 1/9;

sigE = 6.5*rI/sqrt(2)*dx; %0.5; %These are the Gs sigmas; Sqrt(2) is there cause he defines his gaussian's sigmas with the 2. His notation- exp[-(x/sigma)^2]. So multiply by 1/sqrt(2) to equate our notations
sigI = sqrt(1/aI)*sigE; 
sigf = 0.01; % this is the Gf sigma

sEc = (1/sqrt(2))*rc*D_arbor/2*dx;
sIc = (sqrt(1/ac))*sEc; %similar reason for sqrt(2) factor. 

sEc2 = sEc/3;
sIc2 = sIc/3;

%=========================================================================
Inp.L = L;
Inp.d = d;
Inp.nDim = nDim;
Inp.nDimV1 = nDimV1;
Inp.nDimLGN = nDimLGN;
Inp.sigE = sigE;
Inp.sigI = sigI;
Inp.sigf = sigf;
Inp.sEc = sEc;
Inp.sIc = sIc;
Inp.sEc2 = sEc2;
Inp.sIc2 = sIc2;
Inp.ac = ac;
Inp.aE = aE;
Inp.aI = aI;
Inp.af = af;
Inp.a_s = a_s;
Inp.spacing = dx;
Inp.dx = dx; 
Inp.PBC = PBC; %determines if we're using Periodic Boundary conditions (PBC = 1) or not (PBC = 0)
Inp.NORM = NORM;
Inp.ARBOR = ARBOR;
Inp.MILL94 = MILL94;
Inp.nV1 = nV1;
Inp.nThal = nThal;%floor(nV1/2);
Inp.fftpath = fftpath;
Inp.IntMethod = IntMethod;
Inp.c_a = c_a; 

switch IntMethod
    case 'Euler'
        Inp.T = T;
        Inp.Nsteps = round(T/dt);
end
%=========================================================================
%get distances on the grid (in d==1, distV1 gives distances between all pairs of grid-points, in d==2, it gives distances of all grid-points with grid-point at (-L,-L):
switch fftpath
    case 'yes'
        useGdist = 0;
    case 'brute-force'
        useGdist = 1; %useGlobalDistance
end

[distV1sq, Xdist] = neurDistMill(nDim, d, L,PBC, useGdist); %putting a random number as the 5th variable uses globaldist to find the distance. So distV1sq is then nV1 x nThal instead of nDimV1 x nDimLGN
distV1 = sqrt(distV1sq);
Inp.distV1 = distV1;
Inp.Xdist = Xdist; 


[gdist, ~] = neurDistMill(nDim, d, L, PBC, 1); %this way I can use gdist for creating the W so that it's big enough 
gdist = sqrt(gdist);
Inp.gdist = gdist;
%--------------------------------------------------------

if ARBOR 
    [A,~,~,RFinds,Nx_arbor] = Arbor2(Xdist, R_arbor,d, nDim,nDim,L, c_a, 0); %R_arbor = radius of the arbor, c_a = scale factor for the other circle (see Miller 94)
    
    Aon = A;
    Aoff = A;
    
elseif ~ARBOR
    A = ones(nV1,nThal);
    Nx_arbor = nThal;
    RFinds = (1:nV1*nThal)';
end

nRF = sum(Aon(:)~=0) + sum(Aoff(:) ~=0);

%-------------------------------------------------------
Inp.A = A;
Inp.Aon = Aon;
Inp.Aoff = Aoff;
Inp.R_arbor = R_arbor;
Inp.nRF = nRF;
Inp.Nx_arbor = Nx_arbor;
Inp.RFinds = RFinds;
%--------------------------------------------------------
end