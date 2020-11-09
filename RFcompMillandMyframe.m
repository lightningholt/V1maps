function RFcompMillandMyframe(Results, Mill_Results, Inp)

thetas = linspace(0,90, 18);
phis = linspace(0,180, 18);


nDimV1 = Inp.nDimV1;
Nx_arbor = Inp.Nx_arbor;
lw =2 

WRFon = Results.WRFon;
WRFoff = Results.WRFoff;

WRFploton = reshape(WRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFploton = permute(WRFploton, [3,1,4,2]);
WRFploton = reshape(WRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

WRFplotoff = reshape(WRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFplotoff = permute(WRFplotoff, [3,1,4,2]);
WRFplotoff = reshape(WRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

if length(Mill_Results.spatialCorrRF) ~= length(Mill_Results.DThetaNeighbors)
    Mill_Results.spatialCorrRF = RFcorr(Inp, Mill_Results.WRFon, Mill_Results.WRFoff);
end


if length(Results.spatialCorrRF) ~= length(Results.DThetaNeighbors)
    Results.spatialCorrRF = RFcorr(Inp, Results.WRFon, Results.WRFoff);
end

mWRFon = Mill_Results.WRFon;
mWRFoff = Mill_Results.WRFoff;
millWRFploton = reshape(mWRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
millWRFploton = permute(millWRFploton, [3,1,4,2]);
millWRFploton = reshape(millWRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

millWRFplotoff = reshape(mWRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
millWRFplotoff = permute(millWRFplotoff, [3,1,4,2]);
millWRFplotoff = reshape(millWRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);


edges = linspace(-1, 1, 40);

figure(27)

subplot(4,2, [5,7])
histogram(Results.spatialCorrRF, edges, 'Normalization', 'probability');
hold on
histogram(Mill_Results.spatialCorrRF,edges, 'Normalization', 'probability')
hold off
% % axis square
legend('Dynamical Model', 'Prev Model')

x_edge = 14-1;
y_edge = 14-1;
zoom_inds = Nx_arbor:6*Nx_arbor;
midpt = length(zoom_inds/2);

subplot(4,2, 1)
imagesc(WRFploton(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds) - WRFplotoff(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds))
axis square
set(gca, 'Colormap',gray)
%colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
title('Dynamical RFs')

subplot(4,2, 3)
imagesc(millWRFploton(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - millWRFplotoff(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds))
axis square
set(gca, 'Colormap',gray)
%colorbar('Ticks', [min(min(mWRFon-mWRFoff)),max(max(mWRFon-mWRFoff))], 'TickLabels', {'Off', 'On'})
title('Miller RFs')


subplot(4,2,2)
imagesc(Results.zThetaPref)
axis square
hold on
rectangle('Position',[x_edge+zoom_inds(1)/Nx_arbor, y_edge+zoom_inds(1)/Nx_arbor,  zoom_inds(end)/Nx_arbor, zoom_inds(end)/Nx_arbor],...
    'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
hold off
set(gca, 'Colormap', twilight)
colorbar
title('Dynamical Ori Pref')

subplot(4,2,4)
imagesc(Mill_Results.zThetaPref)
axis square
hold on
rectangle('Position',[x_edge+zoom_inds(1)/Nx_arbor, y_edge+zoom_inds(1)/Nx_arbor, zoom_inds(end)/Nx_arbor, zoom_inds(end)/Nx_arbor],...
    'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
hold off
set(gca, 'Colormap',twilight)
title('Miller Ori Pref')
colorbar


subplot(4,2,6)
histogram(Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
legend('Neighbors', 'Random')
title('Dyn Ori Pref')

subplot(4,2,8)
histogram(Mill_Results.DThetaNeighbors,thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(Mill_Results.DThetaRandos,thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
title('Miller Ori Pref')


% =====================



figure(22)
 

subplot(2,4,1)
imagesc(Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Dynamical Ori Pref')

subplot(2,4,2)
histogram(Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
legend('Neighbors', 'Random')

subplot(2,4,3)
imagesc(Mill_Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Miller Ori Pref')

subplot(2,4,4)
histogram(Mill_Results.DThetaNeighbors,thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(Mill_Results.DThetaRandos,thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square


subplot(2,4,5)
imagesc(Results.zPhiPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Dynamical Phase Pref')

subplot(2,4,6)
medPhi = median(Results.DPhiNeighbors);
g = histogram(Results.DPhiNeighbors, phis, 'Normalization', 'probability', 'FaceColor', 'k');
hold on
plot(medPhi, 1.2*max(g.Values))
histogram(Results.DPhiRandos, phis, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square

subplot(2,4,7)
imagesc(Mill_Results.zPhiPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Miller Phase Pref')

subplot(2,4,8)
medPhi = median(Mill_Results.DPhiNeighbors);
g = histogram(Mill_Results.DPhiNeighbors,phis,  'Normalization', 'probability', 'FaceColor', 'k');
hold on
plot(medPhi, 1.2*max(g.Values))
histogram(Mill_Results.DPhiRandos,phis,  'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red',  'LineWidth',lw)
hold off
axis square

f = figure(101);
% f.Units = 'Normalized';
% f.Position =  [0.5, 0.5, 0.7, 0.35];
fontsize=16
subplot(1,2, 2)
medRF = median(Results.spatialCorrRF);
g = histogram(Results.spatialCorrRF, edges,'FaceColor',[0.4660 0.6740 0.1880], 'Normalization', 'probability');
hold on
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility', 'off', 'Color',[0.4660 0.6740 0.1880], 'MarkerSize',12)
medRF = median(Mill_Results.spatialCorrRF);
g = histogram(Mill_Results.spatialCorrRF,edges, 'Normalization', 'probability', 'FaceColor', [0.6350 0.0780 0.1840]);
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility','off', 'Color',[0.6350 0.0780 0.1840], 'MarkerSize',12)
hold off
axis square
legend({'Dynamic V1', 'Static V1'}, 'Box', 'off','Location', 'northwest', 'FontSize', fontsize-2)
ylabel('Probability','FontSize', fontsize)
xlabel('Relative spatial correaltion', 'FontSize', fontsize)
text(-0.2, 1, 'B', 'Units','normalized', 'FontSize', fontsize)
box off;

subplot(1,2,1)
medRF = median(Results.absSpatCorr);
g = histogram(Results.absSpatCorr, edges,'FaceColor',[0.4660 0.6740 0.1880], 'Normalization', 'probability');
hold on
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility', 'off', 'Color',[0.4660 0.6740 0.1880], 'MarkerSize',12)
medRF = median(Mill_Results.absSpatCorr);
g = histogram(Mill_Results.absSpatCorr,edges, 'Normalization', 'probability', 'FaceColor', [0.6350 0.0780 0.1840]);
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility','off', 'Color',[0.6350 0.0780 0.1840], 'MarkerSize',12)
hold off
axis square
legend({'Dynamic V1', 'Static V1'}, 'Box', 'off','Location', 'northwest', 'FontSize', fontsize-2)
ylabel('Probability','FontSize', fontsize)
xlabel('Absolute spatial correaltion', 'FontSize', fontsize)
text(-0.2, 1, 'A', 'Units','normalized', 'FontSize', fontsize)
box off;




