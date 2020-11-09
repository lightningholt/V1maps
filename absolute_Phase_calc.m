% absolute phase comparison
function Results = absolute_Phase_calc(Results, Inp);

nDim = Inp.nDim;
nDimLGN = Inp.nDimLGN;
nDimV1 = Inp.nDimV1;
d = Inp.d;
L = Inp.L;
PBC = Inp.PBC;
Nx_arbor = Inp.Nx_arbor;
R_arbor = Inp.R_arbor;
nV1 = Inp.nV1;
dx = Inp.spacing;

%find the distances
% [distV1sq, Xarbor] = neurDistMill(Nx_arbor, d, R_arbor+0.5,PBC); %R_arbor+0.5 to grab the cell in the middle

Xdist = Inp.Xdist;
Ydist = Inp.Xdist';
% Yarbor = repmat(Xarbor(:,1)',Nx_arbor,1);
k = Nx_arbor/10:Nx_arbor/10:Nx_arbor;
k = 1./k; 
theta = 0:10:170;
theta = theta*pi/180; % to get theta into radians
phi = 0:10:350;
phi = phi*pi/180; 


W = Results.Won - Results.Woff;

RR = zeros(length(phi), length(k), length(theta), nV1);


for th = 1:length(theta)
    for j = 1:length(k)
        for p = 1:length(phi)
            k_x = k(j)*cos(theta(th)); %finds the projection of k in the x-direction
            k_y = k(j)*sin(theta(th));
            grating = sin((k_x*Xdist + k_y*Ydist)*2*pi + phi(p));
            RR(p,j,th,:) = W*grating(:);
        end
    end
end


mTheta = zeros(nV1,1);
mK = mTheta;
mPhi = mTheta;


for nn = 1:nV1
    Resp = RR(:,:,:,nn);
    [~, idx] = max(Resp(:));
    [mPhi(nn), mK(nn), mTheta(nn)] = ind2sub([length(phi), length(k), length(theta)], idx);
end

PrefPhi = phi(mPhi)*180/pi;
PrefK = k(mK);
PrefTheta = theta(mTheta)*180/pi;


PrefK = reshape(PrefK, nDimV1, nDimV1);
PrefTheta = reshape(PrefTheta, nDimV1, nDimV1);
PrefPhi = reshape(PrefPhi, nDimV1, nDimV1);

for nn = 1:nV1
%     rTheta(:, nn) = RR(mPhi(nn), mK(nn), :, nn); %response at max phase and spatial frequency (theta is the only variable remaining)
    rPhi (:, nn) = RR(:, mK(nn), mTheta(nn), nn); % response at max theta and spatial frequency (phi si the only variable remaining)
end

zPhi = exp(1i*phi)*rPhi;%./(sum(abs(rPhi))); %Miller normalizes by rms(abs(rPhi))*length(phi)
zPhiSel = abs(zPhi)./(rms(rPhi)*length(phi));
zPhiPref = angle(zPhi);

zAbsPhiSel = reshape(zPhiSel, nDimV1, nDimV1);
zAbsPhiPref = reshape(zPhiPref, nDimV1, nDimV1)*180/pi;
zAbsPhiPref(zPhiPref < 0) = 360 + zPhiPref(zPhiPref <0);


%find difference between neighbors and randos
X = zeros(nDimV1^d, 1); %for periodic boundary conditions, these are bounded between 1 and L/2
Y = zeros(nDimV1^d, 1);


for i = 1:nDimV1^d
    [X(i), Y(i)] = ind2sub([nDimV1, nDimV1], i);
end

XX = X;
YY = Y;

X = X-1;
Y = Y-1;
R = sqrt(2)+0.1;

DAbsPhiNeighbors = [];
DAbsPhiRandos = [];

for i = 1:nV1
    for j = (i+1):nV1
        if Inp.PBC
            if abs(X(i) - X(j)) > Inp.nDim/2
                xdist = Inp.nDim - abs(X(i) - X(j));
            else
                xdist = X(i) - X(j);
            end
            
            if abs(Y(i) - Y(j)) > Inp.nDim/2
                ydist = Inp.nDim - abs(Y(i) - Y(j));
            else
                ydist = Y(i) - Y(j);
            end
        else
            xdist = X(i) - X(j);
            ydist = Y(i) - Y(j);
            
        end
        
        if xdist^2 + ydist^2 <= R^2 %need to use the PBC X and Y here.
%         if (X(i) - X(j))^2 + (Y(i) - Y(j))^2 <= R^2 %need to use the PBC X and Y here.

%             DThetaNeighbors = [DThetaNeighbors; abs(PrefTheta(XX(i), YY(i)) - PrefTheta(XX(j), YY(j)))];
            
            DAbsPhiNeighbors = [DAbsPhiNeighbors; abs(zAbsPhiPref(XX(i), YY(i)) - zAbsPhiPref(XX(j), YY(j)))];
            
            randX = randi(nDimV1 , 1);
            randY = randi(nDimV1, 1); 
            
            
%             DThetaRandos = [DThetaRandos; abs(PrefTheta(XX(i), YY(i)) - PrefTheta(randX, randY))];
            DAbsPhiRandos = [DAbsPhiRandos; abs(zAbsPhiPref(XX(i), YY(i)) - zAbsPhiPref(randX, randY))];
        end
    end
end
    
Results.zAbsPhiSel = zAbsPhiSel;
Results.zAbsPhiPref = zAbsPhiPref;
Results.DAbsPhiNeighbors = DAbsPhiNeighbors;
Results.DAbsPhiRandos = DAbsPhiRandos;

end