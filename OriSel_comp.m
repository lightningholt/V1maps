% compare the orientation preference between neighboring cells for
% different models

function OriSel_comp(Big_Results, Inp, varargin)

if nargin > 2
    cmax = varargin{1};
else
    cmax = max(max(Big_Results.Res.zThetaSel));
    cmax1 = max(max(Big_Results.Mill_Res.zThetaSel));
    if cmax < cmax1 
        cmax = cmax1;
    end 
end

cmin = 0;

% Big_Results should contain the results for the different conditions, i.e.
% Res - the cortical dynamic results, Mill_Res - the results using Miller
% 94's formalism, the different response kernel Res, ell- using elliptical
% rather than circular MH's for the K's, rot_ - each cortical cell has a
% randomly rotated ellipse, squi_ - each cortical cell has a randomly
% squished ellipse, squi_rot- each cortical has a randomly squished and
% rotated ellipse. 

%Inp the struct containing all the inputs
% ii - the filenumber 
% small_comp - just compare Mill_Res and Res -- the ellispses are all
% basically the same and don't look any bette than Res.


Results = Big_Results.Res;
Mill_Results =Big_Results.Mill_Res ;
% ell_Results = Big_Results.ell_Res ;
% rot_Results = Big_Results.rot_Res;
% squi_Results = Big_Results.squi_Res ;
% squi_rot_Results = Big_Results.squi_rot_Res;

[millWRFploton, millWRFplotoff] = readyWRFtoPlot(Mill_Results, Inp);
[WRFploton, WRFplotoff] = readyWRFtoPlot(Results, Inp);

% figname = ['CompOri_',num2str(ii)];

figure(999)
thetas = 0:5:90;
lw = 2.5;
%{

subplot(2,3,1)
imagesc(Mill_Results.zThetaSel)
axis square
set(gca, 'Colormap', redblue)
colorbar
caxis([cmin cmax])
title('Miller')


subplot(2,3,2)
imagesc(Results.zThetaSel);
axis square
set(gca, 'Colormap', redblue)
colorbar
caxis([cmin, cmax])
title('Dyn')

subplot(2,3,3)
imagesc(ell_Results.zThetaSel);
axis square
set(gca, 'Colormap', redblue)
colorbar
caxis([cmin, cmax])
title('Ell+Dyn')


subplot(2,3,4)
imagesc(rot_Results.zThetaSel);
axis square
set(gca, 'Colormap', redblue)
colorbar
caxis([cmin, cmax])
title('rot+Dyn')



subplot(2,3,5)
imagesc(squi_Results.zThetaSel);
axis square
set(gca, 'Colormap', redblue)
colorbar
caxis([cmin, cmax])
title('squi+Dyn')



subplot(2,3,6)
imagesc(squi_rot_Results.zThetaSel);
axis square
set(gca, 'Colormap', redblue)
colorbar
caxis([cmin, cmax])
title('squi+rot+Dyn')
%}

x_edge = 14-1;
y_edge = 14-1;
Nx_arbor = Inp.Nx_arbor;
zoom_inds = Inp.Nx_arbor:6*Inp.Nx_arbor;
midpt = length(zoom_inds/2);

figure(1000)
thetas = 0:5:90;
lw = 2.5;
fontsize= 14;

c1 = max(Results.zThetaSel(:));
c2 = max(Mill_Results.zThetaSel(:));

if c1 > c2
    cmax = c1;
else
    cmax = c2;
end


subplot(2,2, 1)
imagesc(WRFploton(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds) - WRFplotoff(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds))
axis square
set(gca, 'Colormap',gray)
text(-0.15, 1, 'A', 'Units','normalized', 'FontSize', fontsize)
%colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
%title('Dynamical RFs')
axis off

subplot(2,2, 3)
imagesc(millWRFploton(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - millWRFplotoff(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds))
axis square
set(gca, 'Colormap',gray)
text(-0.15, 1, 'C', 'Units','normalized', 'FontSize', fontsize)
%colorbar('Ticks', [min(min(mWRFon-mWRFoff)),max(max(mWRFon-mWRFoff))], 'TickLabels', {'Off', 'On'})
%title('Miller RFs')
axis off



subplot(2,2,2)
%imagesc(Results.zThetaSel);
% set(gca, 'Colormap', redblue)
% c = colorbar;
% c.Label.String = 'Orientation selectivity';
% c.Label.FontSize = 14
% caxis([cmin, cmax])
% xticks([])
% yticks([])
% hold on
% rectangle('Position',[x_edge+zoom_inds(1)/Nx_arbor, y_edge+zoom_inds(1)/Nx_arbor,  zoom_inds(end)/Nx_arbor, zoom_inds(end)/Nx_arbor],...
%     'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
% hold off
histogram(Results.zThetaSel, 'FaceColor',[0.4660 0.6740 0.1880], 'Normalization', 'probability')
axis square
median(Results.zThetaSel(:))
meanSel = mean(Results.zThetaSel(:))
text(-0.3, 1, 'B', 'Units','normalized', 'FontSize', fontsize)
hold on
xline(meanSel, 'LineWidth', 3, 'LineStyle', '--', 'Color', [0.4660 0.6740 0.1880])
hold off
box off 

subplot(2,2,4)
% imagesc(Mill_Results.zThetaSel)
% set(gca, 'Colormap', redblue)
% c = colorbar;
% c.Label.String = 'Orientation selectivity';
% c.Label.FontSize = 14
% caxis([cmin cmax])
% box on
% xticks([])
% yticks([])
% hold on
% rectangle('Position',[x_edge+zoom_inds(1)/Nx_arbor, y_edge+zoom_inds(1)/Nx_arbor,  zoom_inds(end)/Nx_arbor, zoom_inds(end)/Nx_arbor],...
%     'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
% hold off
histogram(Mill_Results.zThetaSel, 'Normalization', 'probability', 'FaceColor', [0.6350 0.0780 0.1840]);
axis square
xlabel('Orientation selectivity index');
median(Mill_Results.zThetaSel(:))
meanSel = mean(Mill_Results.zThetaSel(:))
text(-0.3, 1, 'D', 'Units','normalized', 'FontSize', fontsize)
hold on
xline(meanSel,  'LineWidth', 3, 'LineStyle', '--', 'Color', [0.6350 0.0780 0.1840])
hold off
box off
% title('Miller')




end