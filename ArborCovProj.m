%Function to find the covariance (correlation) within the arbor

function [wp, CL] = ArborCovProj(Inp, Won, Woff)

nThal = Inp.nThal; %
R = Inp.R_arbor; %radius of the arbor
nDim = Inp.nDimV1; 
dx = 2*Inp.L/nDim;
Nrf = length(-R:dx:R); %Number/dim in RF at the max (across the middle)
Indrf = ceil(-Nrf/2):floor(Nrf/2);

d = Inp.d;
% Neff = Nrf^d;

%calculate the distance between cells within the arbor
[xs0, ys0] = meshgrid(-R:dx:R);
% Nrf = size(xs0,1);

%now find the distances that correspond with the cells in the arbor
%to do that initialize the arbor
A = Inp.A;
Amask = reshape(A(1,:)>0, nDim, nDim);
Amask2 = A > 0;

Cnn = Inp.Cnn;
Cnn_fast = Inp.Cnn_fast;

%{
%fftshift will put the first neurons arbor in the middle of the grid
A_1 = fftshift(reshape(Amask2(1,:), nDim, nDim));
A_1 = A_1(nDim/2+Indrf+1, nDim/2+Indrf+1);
%double check

figure(54)
imagesc(A_1);

% zero2D = xs0(:) * 0;

Cnn2 = fftshift(Cnn);
Cnn2 = Cnn2(nDim/2+Indrf+1, nDim/2+Indrf+1);
Cnn2 = Cnn2(A_1);

Cnn_fast2 = fftshift(Cnn_fast);
Cnn_fast2 = Cnn_fast2(nDim/2+Indrf+1, nDim/2+Indrf+1);
Cnn_fast2 = Cnn_fast2(A_1);

dC_slow = 1.5*Cnn2;
dC_fast = 1.5*Cnn_fast2;

[~, v1] = svds(dC_slow, 1);
[~, v2] = svds(dC_fast, 1);

if v1 ~= 1
    Cnn = (1/v1) * Cnn;
    Cnn2 = (1/v1) * Cnn2;
    dC_slow = 1.5*Cnn2;
    [~, v1] = svds(dC_slow,1);
    v1
end

if v2 ~= 1
    Cnn_fast = (1/v2) * Cnn_fast;
    Cnn_fast2 = (1/v2) * Cnn_fast2;
    dC_fast = 1.5*Cnn_fast2;
    [~, v2] = svds(dC_fast,1);
    v2
end

%}

e1 = Inp.e1RF;
e2 = Inp.e2RF;
e3 = Inp.e3RF;


if e2'*e1 > 0.01
    warning('e1 and e2 are not orthogonal yo')
    e2'*e1
elseif e2'*e3 > 0.01
    warning('e2 and e3 could be more orthogonal')
    e2'*e3
elseif e1'*e3 > 0.01
    warning('e1 and e3 could be more orthogonal')
    e1'*e3
end


figure(55)
plot(1:length(e1), e1, 1:length(e1), e2)
legend('Slow Feature', 'fast Feature')

Neff = sum(Amask2(1,:));

W = zeros(Inp.nV1, Neff);

for ii = 1:Inp.nV1
    W(ii,:) = circshift(Won(ii,Amask2(ii,:)) - Woff(ii, Amask2(ii,:)), ii -1);
end

wp1 = W*e1;
slow_mean = mean(wp1.^2);
wp1 = reshape(wp1, Inp.nDimV1, Inp.nDimV1);
wp1_roughness = nearNeighbor2D(wp1, Inp)/slow_mean;
CL1 = correlationLengthCalc(wp1, Inp);

wp2 = W*e2;
fast_mean = mean(wp2.^2);
wp2 = reshape(wp2, Inp.nDimV1, Inp.nDimV1);
wp2_roughness = nearNeighbor2D(wp2, Inp)/fast_mean;
CL2 = correlationLengthCalc(wp2, Inp);

wp3 = W*e3;
neither_mean = mean(wp3.^2);
wp3 = reshape(wp3, Inp.nDimV1, Inp.nDimV1);
wp3_roughness = nearNeighbor2D(wp3, Inp)/neither_mean;
CL3 = correlationLengthCalc(wp3, Inp);

expectedScaledDifference = (max(Inp.lambda1*abs(e1)))/(max(Inp.lambda2*abs(e2)))
actualScaledDifference = max(abs(wp1(:)))/max(abs(wp2(:)))

n_space = 1:length(CL1);

% figure(67);
% subplot(1,3,1)
% imagesc(wp1)
% title("Selectivity to slow C's")
% botstr = strcat('Mean Selectivity =', num2str(slow_mean));
% botstr = [botstr newline ' Roughness = ' num2str(wp1_roughness)];
% xlabel(botstr);
% % caxis([cmin, cmax])
% colorbar
% 
% subplot(1,3,2)
% imagesc(wp2)
% title("Selectivity to fast C's")
% botstr = strcat('Mean Selectivity =', num2str(fast_mean));
% botstr = [botstr newline ' Roughness = ' num2str(wp2_roughness)];
% xlabel(botstr);
% % caxis([cmin, cmax])
% colorbar
% 
% subplot(1,3,3)
% imagesc(wp3)
% title("Selectivity to neither")
% botstr = strcat('Mean Selectivity =', num2str(neither_mean));
% botstr = [botstr newline ' Roughness = ' num2str(wp3_roughness)];
% xlabel(botstr);
% % caxis([cmin, cmax])
% colorbar
% 
% figure(97)
% %exportpdf(97, 'LinearC_organizations_46', [8,0.5])
% subplot(2,3,1)
% imagesc(wp1)
% title("Selectivity to slow C's")
% botstr = strcat('Mean Selectivity =', num2str(slow_mean));
% botstr = [botstr newline ' Roughness = ' num2str(wp1_roughness)];
% xlabel(botstr);
% % caxis([cmin, cmax])
% colorbar
% 
% subplot(2,3,2)
% imagesc(wp2)
% title("Selectivity to fast C's")
% botstr = strcat('Mean Selectivity =', num2str(fast_mean));
% botstr = [botstr newline ' Roughness = ' num2str(wp2_roughness)];
% xlabel(botstr);
% % caxis([cmin, cmax])
% colorbar
% 
% subplot(2,3,3)
% imagesc(wp3)
% title("Selectivity to neither")
% botstr = strcat('Mean Selectivity =', num2str(neither_mean));
% botstr = [botstr newline ' Roughness = ' num2str(wp3_roughness)];
% xlabel(botstr);
% % caxis([cmin, cmax])
% colorbar
% 
% subplot(2,3,4)
% fftCslow = abs(fft2(dC_slow));
% imagesc(fftCslow)
% title('fft2(Cslow)')
% colorbar
% 
% subplot(2,3,5)
% fftCfast = abs(fft2(dC_fast));
% imagesc(fftCfast)
% title('fft2(Cfast)')
% colorbar;
% 
% subplot(2,3, 6)
% plot(n_space, [CL1, CL2, CL3], "LineWidth", 2.25)
% legend('Slow feature', 'Fast Feature', 'Neither')


wp(:,:, 1) = wp1;
wp(:,:,2) = wp2;
wp(:,:,3) = wp3;

CL = [CL1, CL2, CL3];


    
