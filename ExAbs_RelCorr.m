%to make a cartoon example of absolute and spatial Correlations
% and show historgrams

Results = Big_Results.Res;
Mill_Results = Big_Results.Mill_Res;

fontsize=12
W = Big_Results.Res.Won - Big_Results.Res.Woff;
% W1 = fftshift(reshape(W(531,:), Inp.nDimLGN, Inp.nDimLGN));
% W2 = fftshift(reshape(W(563,:), Inp.nDimLGN, Inp.nDimLGN));

W1 = reshape(W(531,:), Inp.nDimLGN, Inp.nDimLGN);
W2 = reshape(W(563,:), Inp.nDimLGN, Inp.nDimLGN);


Wrf = Big_Results.Res.WRFon - Big_Results.Res.WRFoff;
Wrf1 = (reshape(Wrf(1,:), Inp.Nx_arbor, Inp.Nx_arbor));
Wrf2 = (reshape(Wrf(2,:), Inp.Nx_arbor, Inp.Nx_arbor));


figure(89)
subplot(3,4,1)
imagesc(W1)
colormap('gray')
axis square
xticks([])
yticks([])
%text(-0.2, 1, 'A', 'Units','normalized', 'FontSize', fontsize)

subplot(3,4,2)
imagesc(W2)
colormap('gray')
axis square
xticks([])
yticks([])

subplot(3,4,3)
imagesc(Wrf1)
colormap('gray')
axis square
xticks([])
yticks([])
%text(-0.2, 1, 'B', 'Units','normalized', 'FontSize', fontsize)

subplot(3,4,4)
imagesc(Wrf2)
colormap('gray')
axis square
xticks([])
yticks([])

edges = -1:0.05:1;

subplot(3,4, [5,6, 9,10])%, 13,14])
medRF = median(Results.absSpatCorr);
g = histogram(Results.absSpatCorr, edges,'FaceColor',[0.4660 0.6740 0.1880], 'Normalization', 'probability');
% [ResAbs, Edges] = histcounts(Results.absSpatCorr, edges, 'Normalization', 'probability');
% g = barh(Edges(2:end), ResAbs, 'FaceColor', [0.4660 0.6740 0.1880]);
hold on
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility', 'off', 'Color',[0.4660 0.6740 0.1880], 'MarkerSize',12)
medRF = median(Mill_Results.absSpatCorr);
g = histogram(Mill_Results.absSpatCorr,edges, 'Normalization', 'probability', 'FaceColor', [0.6350 0.0780 0.1840]);
% [MillAbs, Edges] = histcounts(Mill_Results.absSpatCorr,edges, 'Normalization', 'probability');
% g = barh(Edges(2:end), MillAbs, 'FaceColor', [0.6350 0.0780 0.1840]);
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility','off', 'Color',[0.6350 0.0780 0.1840], 'MarkerSize',12)
hold off
% axis square
legend({'Dynamic V1', 'Static V1'}, 'Box', 'off','Location', 'northwest', 'FontSize', fontsize-2)
%text(-0.2, 1, 'C', 'Units','normalized', 'FontSize', fontsize)
box off;
% set(gca, 'view', [90, -90]);
ylabel('Probability','FontSize', fontsize)
xlabel('Absolute spatial correaltion', 'FontSize', fontsize)

subplot(3,4,[7,8,11,12])%,15,16])
medRF = median(Results.spatialCorrRF);
g = histogram(Results.spatialCorrRF, edges,'FaceColor',[0.4660 0.6740 0.1880], 'Normalization', 'probability');
hold on
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility', 'off', 'Color',[0.4660 0.6740 0.1880], 'MarkerSize',12)
medRF = median(Mill_Results.spatialCorrRF);
g = histogram(Mill_Results.spatialCorrRF,edges, 'Normalization', 'probability', 'FaceColor', [0.6350 0.0780 0.1840]);
plot(medRF, 1.1*max(g.Values), 'v', 'HandleVisibility','off', 'Color',[0.6350 0.0780 0.1840], 'MarkerSize',12)
hold off
% axis square
legend({'Dynamic V1', 'Static V1'}, 'Box', 'off','Location', 'northwest', 'FontSize', fontsize-2)
ylabel('Probability','FontSize', fontsize)
xlabel('Relative spatial correaltion', 'FontSize', fontsize)
%text(-0.2, 1, 'D', 'Units','normalized', 'FontSize', fontsize)
box off;

corrcoef(W1, W2)
corrcoef(Wrf1, Wrf2)
