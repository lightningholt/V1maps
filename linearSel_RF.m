% make linearC plots + RF plot in the same figure

%IMPORTANTLY: This is assuming that Cnn = Cff and Cnf = Cfn
function Results = linearSel_RF(Won, Woff, WRFon, WRFoff, Inp, Results)

[WRFploton, WRFplotoff] = RF_plots(WRFon, WRFoff, Inp);

if and(isfield(Results, "wp"), isfield(Results, "CL"))
    wp = Results.wp;
    CL = Results.CL;
else
    [wp, CL] = linear2Dselectivity(Inp, Won, Woff);
    Results.wp = wp;
    Results.CL = CL;
end

if and(isfield(Results, "wpRF"), isfield(Results, "CLRF"))
    wpRF = Results.wpRF;
    CLRF = Results.CLRF;
    e1RF = Results.e1RF;
    e2RF = Results.e2RF;
    e3RF = Results.e3RF;
else
    [wpRF, CLRF, e1RF, e2RF, e3RF] = ArborCovProj(Inp, Won, Woff);
    Results.wpRF = wpRF;
    Results.CLRF = CLRF;
    Results.e1RF = e1RF;
    Results.e2RF = e2RF;
    Results.e3RF = e3RF;
end

% wp/wpRF are both (neurons, neurons, x) where x is in [1, 2, 3] and corresponds to slow, fast, neither
% Similarly CL/CLRF are [neural space, x]

e1 = Inp.e1;
e2 = Inp.e3;
e3 = Inp.e3;

% In order to make e1RF have dimensions that I can display next to the
% others
nV1 = Inp.nV1; %
R = Inp.R_arbor; %radius of the arbor
nDim = Inp.nDimV1; 
dx = 2*Inp.L/nDim;
Nrf = length(-R:dx:R);
Indrf = ceil(-Nrf/2):floor(Nrf/2);

d = Inp.d;
% Neff = Nrf^d;

%calculate the distance between cells within the arbor
[xs0, ys0] = meshgrid(-R:dx:R);
% Nrf = size(xs0,1);

%now find the distances that correspond with the cells in the arbor
%to do that initialize the arbor
A = Inp.A;
Amask = A>0;

%fftshift will put the first neurons arbor in the middle of the grid
A_1 = fftshift(reshape(Amask(1,:), nDim, nDim));
A_1 = A_1(nDim/2+Indrf+1, nDim/2+Indrf+1);

slow_feature = zeros(Nrf, Nrf);
fast_feature = zeros(Nrf, Nrf);
arbor_inds = find(A_1);

for ii = 1:length(arbor_inds)
    slow_feature(arbor_inds(ii)) = e1RF(ii);
    fast_feature(arbor_inds(ii)) = e2RF(ii);
end
    
n_space = (1:length(CL))*Inp.dx;

figure(17);

ax1 = subplot(4, 3, [1:2, 4:5]);
imagesc(WRFploton - WRFplotoff)
colormap(ax1, 'gray');
colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
title('Difference between On & Off')

ax2 = subplot(4,3,3);
imagesc(wpRF(:,:,1)) % C_slow proj
colormap(ax2, 'cool')
colorbar;
title('C_{slow} Mapping')

ax3 = subplot(4,3,6);
% imagesc(reshape(e1, Inp.nDimV1, Inp.nDimV1))
imagesc(slow_feature)
title('C_{slow} Feature')
colormap(ax3, 'gray')
colorbar;

ax4 = subplot(4,3, 7);
% plot(n_space, CL1, n_space, CL2, n_space, CL3, 'LineWidth', 2.5);
plot(n_space, CLRF,'LineWidth', 2.5);
legend('slow feature', 'fast feature', 'neither', 'Slow Corr','Fast Corr')
title('Correlation lengths')
xlabel('Neural distance')

ax5 = subplot(4,3,9);
imagesc(wpRF(:,:,2))
colormap(ax5, 'cool')
colorbar;
title('C_{fast} Mapping')

ax6 = subplot(4,3,8);
% imagesc(reshape(e2, Inp.nDimV1, Inp.nDimV1))
imagesc(fast_feature)
colormap(ax6, 'gray')
title('C_{fast} Feature')
colorbar;

RF_cell1 = WRFploton(1:Inp.Nx_arbor, 1:Inp.Nx_arbor) - WRFplotoff(1:Inp.Nx_arbor, 1:Inp.Nx_arbor);
% Lin_cell1 = wp1(1).*reshape(e1, Inp.nDimV1, Inp.nDimV1) + wp2(1).*reshape(e2, Inp.nDimV1, Inp.nDimV1);
Lin_cell1 = wpRF(1,1,1).*slow_feature + wpRF(1,1,2).*fast_feature * max(abs(e1RF))/max(abs(e2RF));

% figure(18)
ax7 = subplot(4,2,7);
imagesc(RF_cell1)
colorbar;
colormap(ax7, 'gray')
title('RF of Cell 1')

subplot(4,2,8)
imagesc(Lin_cell1)
colorbar;
title('Lin Combination of C_{slow} + C_{fast} feature')

figure(18)
ax_rf = subplot(4,4, [1:2, 5:6]);
imagesc(WRFploton - WRFplotoff)
colormap(ax_rf, 'gray');
colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
title('Difference between On & Off')

ax_slow = subplot(6,4,[3,7]);
imagesc(wpRF(:,:,1))
colorbar;
colormap(ax_slow, 'hot')
title('Slow Feature Mapping')
axis equal
axis off

ax_fast = subplot(6,4,[4,8]);
imagesc(wpRF(:,:,2))
colorbar;
colormap(ax_fast, 'hot')
title('Fast Feature Mapping')
axis equal
axis off

ax_cl = subplot(6,4, [11,12]);
plot(n_space, CLRF, 'LineWidth', 2.25);
legend('slow feature', 'fast feature', 'neither', 'Slow Corr','Fast Corr')
title('Correlation lengths')
xlabel('Neural distance')

ax_spat = subplot(4,4, [9:10, 13:14]);
histogram(Results.spatialCorrRF, 'HandleVisibility','off', 'Normalization','probability')
hold on
spat_med = median(Results.spatialCorrRF);
spat_mean = mean(Results.spatialCorrRF);
spat_mean = num2str(spat_mean);
spat_title = ['Correlation Between Neighbors, mean =',spat_mean];
yy = ylim;
plot(spat_med, 0.95*yy(end), 'v', 'LineWidth', 2, 'MarkerSize',15);
% line([spat_mean, spat_mean], yy, 'LineWidth', 2);
legend('Median')
title(spat_title)
hold off

ax_ori = subplot(6,4, [15,19]);
imagesc(Results.zThetaPref)
colorbar;
colormap(ax_ori, 'hsv')
title('Prefered Ori')
axis equal
axis off

ax_phase = subplot(6,4,[16,20]);
imagesc(Results.zPhiPref)
colorbar;
colormap(ax_phase, 'hsv')
title('Prefered Phase')
axis equal
axis off

edges = 0:10:90;

ax_dtheta = subplot(6,4,23);
histogram(Results.DThetaNeighbors,edges,  'Normalization', 'probability')
hold on
histogram(Results.DThetaRandos, edges, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2)
legend('Neighbors', 'Random', 'Location', 'southeast');
title('Ori Difference')
% xlabel('Diff of Preferred Phase')
ylabel('Probability')
hold off

edgesP = 0:10:180;
ax_dphi = subplot(6,4,24);
histogram(Results.DPhiNeighbors,edgesP,  'Normalization', 'probability')
hold on
histogram(Results.DPhiRandos, edgesP, 'Normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2)
legend('Neighbors', 'Random', 'Location', 'southeast');
title('Phase Difference')
% xlabel('Diff of Preferred Phase')
ylabel('Probability')
hold off


