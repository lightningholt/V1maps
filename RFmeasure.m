%Created by Caleb Holt
% to measure the receptive fields of the a 2D grid of on and off cells. 
% multiply the RF with a sine wave that spans the whole 

%Define network parameters
function [PrefTheta, PrefK, PrefPhi, zThetaSel, zThetaPref, zPhiSel, zPhiPref] = RFmeasure(Inp, WRFon, WRFoff)
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

Xarbor = -R_arbor:dx:R_arbor;
Yarbor = Xarbor';
% Yarbor = repmat(Xarbor(:,1)',Nx_arbor,1);
k = Nx_arbor/10:Nx_arbor/10:Nx_arbor;
k = 1./k; 
theta = 0:10:170;
theta = theta*pi/180; % to get theta into radians
phi = 0:10:350;
phi = phi*pi/180; 

%load previous results
onRF = reshape(WRFon',Nx_arbor, Nx_arbor, nDimV1^d); %first two dimensions give on RF of cell at 3rd dimension
offRF = reshape(WRFoff',Nx_arbor, Nx_arbor, nDimV1^d); 
RF = onRF - offRF;

RR = zeros(length(phi), length(k), length(theta), nV1);

for nn = 1:nV1
    for th = 1:length(theta)
        for j = 1:length(k)
            for p = 1:length(phi)
                k_x = k(j)*cos(theta(th)); %finds the projection of k in the x-direction
                k_y = k(j)*sin(theta(th));
                J = RF(:,:,nn).*sin((k_x*Xarbor + k_y*Yarbor)*2*pi + phi(p));
                J = J(:);
                RR(p,j,th,nn) = sum(J);
            end
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

WRFploton = reshape(WRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFploton = permute(WRFploton, [3,1,4,2]);
WRFploton = reshape(WRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

WRFplotoff = reshape(WRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFplotoff = permute(WRFplotoff, [3,1,4,2]);
WRFplotoff = reshape(WRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

PrefK = reshape(PrefK, nDimV1, nDimV1);
PrefTheta = reshape(PrefTheta, nDimV1, nDimV1);
PrefPhi = reshape(PrefPhi, nDimV1, nDimV1);

halfjet = length(jet)/2;
jet_wrap = vertcat(jet, flipud(jet)); %vertcat(jet((1:halfjet),:), flipud(jet((halfjet + 1):end,:));%

% figure(101)
% imagesc(PrefTheta);
% colormap(jet_wrap)
% colorbar
% title('Ori Map')
% 
% figure(102)
% imagesc(PrefK)
% colorbar
% title('Spatial freq Map')
% 
% figure(103)
% imagesc(PrefPhi)
% colormap(jet_wrap)
% colorbar
% title('Phase Preference Map') 

% To analyze the orientation selectivity of the network
for nn = 1:nV1
    rTheta(:, nn) = RR(mPhi(nn), mK(nn), :, nn); %response at max phase and spatial frequency (theta is the only variable remaining)
    rPhi (:, nn) = RR(:, mK(nn), mTheta(nn), nn); % response at max theta and spatial frequency (phi si the only variable remaining)
end

zTheta = exp(2*1i*theta)*rTheta;%./(sum(abs(rTheta))); %Miller divides by rms(abs(rTheta))*length(theta)
zThetaSel = abs(zTheta)./(rms(rTheta)*length(theta));
zThetaPref = 1/2*angle(zTheta);

zPhi = exp(1i*phi)*rPhi;%./(sum(abs(rPhi))); %Miller normalizes by rms(abs(rPhi))*length(phi)
zPhiSel = abs(zPhi)./(rms(rPhi)*length(phi));
zPhiPref = angle(zPhi);

zThetaSel = reshape(zThetaSel, nDimV1, nDimV1);
zThetaPref = reshape(zThetaPref, nDimV1, nDimV1)*180/pi;
zThetaPref(zThetaPref < 0) = 180 + zThetaPref(zThetaPref < 0);
zPhiSel = reshape(zPhiSel, nDimV1, nDimV1);
zPhiPref = reshape(zPhiPref, nDimV1, nDimV1)*180/pi;
zPhiPref(zPhiPref < 0) = 360 + zPhiPref(zPhiPref <0);

figure(104)
imagesc(zThetaSel)
colormap(jet)
colorbar
title('Ori Selectivity')
% 
% figure(105)
% imagesc(zThetaPref)
% colormap(jet_wrap)
% colorbar
% title('Ori Pref 2')
% 
% figure(106)
% imagesc(zPhiSel)
% colormap(jet)
% colorbar
% title('Phase Selectivity')
% 
% figure(107)
% imagesc(zPhiPref)
% colormap(jet_wrap)
% colorbar
% title('Phase Pref 2')

figure(111)
% ax1 = subplot(4, 1, 1);
% imagesc(WRFploton - WRFplotoff)
% colormap(ax1, 'gray')
% title('On/Off Dominance')
% colorbar

ax2 = subplot(3,2,1);
imagesc(PrefTheta);
colormap(ax2, jet)
colorbar
title('Pref Ori Map (Max Resp)')

ax3 = subplot(3,2,2);
imagesc(PrefPhi)
colormap(ax3, jet)
colorbar
title('Pref Phase Map (Max Resp)')

ax4 = subplot(3,2,3);
imagesc(zThetaPref)
colormap(ax4, jet)
colorbar
title('Pref Ori Map - arg(z)')

ax5 = subplot(3,2,4);
imagesc(zPhiPref)
colormap(ax5, jet)
colorbar
title('Pref Phase Map - arg(z)')

ax6 = subplot(3,2,5);
imagesc(zThetaSel)
colormap(ax6, jet)
colorbar
title('Ori Sel Map - abs(z)')

ax7 = subplot(3,2,6);
imagesc(zPhiSel)
colormap(ax7, jet)
colorbar
title('Phase Sel Map - abs(z)')

end
% zTheta = squeeze(RR(mPhi, mK, :,:));
