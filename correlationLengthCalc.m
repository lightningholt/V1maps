%function to find the correlation length of the 2D maps

function CL = correlationLengthCalc(proj, Inp)

%extract the relevant size parameters from the Inp structure
d = Inp.d;
nDimV1 = Inp.nDimV1;
nV1 = Inp.nDimV1;
L = Inp.L; 
Xdist = Inp.Xdist;
c_a = 1; %complete overlap between dendrites and axons (or makes a pill-box arbor);

%find lengths I want to look over
dx = 2*L/nDimV1;
if d == 1
    max_x = L;%nV1*dx;
else
    max_x = 10*dx; % floor(sqrt(2)*L);
end
x = dx:dx:max_x;

%Correlation length initialize
CL = zeros(length(x),1);

for rr = 1:length(x)
    %make an annulus of width 1 neuron by finding arbor at 2 distances
    Arbor_rad = x(rr);
    [Ar,~,~,~,~] = Arbor2(Xdist, Arbor_rad,d, nDimV1,nDimV1,L, c_a, 1);%1 at the end ensures pill-box
    
    if Arbor_rad > dx
        Smaller_rad = x(rr-1);
        [Ar_1,~,~,~,~] = Arbor2(Xdist, Smaller_rad, d, nDimV1, nDimV1, L, c_a,1); %1 at the end ensures pill-box 
    elseif Arbor_rad == dx
        Ar_1 = eye(size(Ar));
    else
        Ar_1 = zeros(size(Ar));
    end  
    
    Ar = logical(Ar - Ar_1); %so now arbor is just annulus 1 neuron thick 
    %find the number of cells in that annulus. 
    nRF = sum(Ar(1,:));
    
    xc = mean(proj(:).*(Ar*proj(:))/nRF); %second term effectively takes the mean of the "RF"
    xs = mean(proj(:)).*mean(Ar*proj(:)/nRF);
    CL(rr) = xc-xs; 
end

 %find the zero displacement correlation, nRF =1;
xc = mean(proj(:).*(eye(size(Ar))*proj(:)/1));
xs = mean(proj(:)).*mean(eye(size(Ar))*proj(:)/1);
CL = [xc-xs; CL]./(xc-xs);
