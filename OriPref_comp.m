% compare the orientation preference between neighboring cells for
% different models

function OriPref_comp(Big_Results, Inp, ii, small_comp)

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

if nargin < 3
    ii = 0
    small_comp = 0
elseif nargin > 3
    small_comp = 1;
end


Results = Big_Results.Res;
Mill_Results =Big_Results.Mill_Res ;

% ell_Results = Big_Results.ell_Res ;
% rot_Results = Big_Results.rot_Res;
% squi_Results = Big_Results.squi_Res ;
% squi_rot_Results = Big_Results.squi_rot_Res;
% 
% figname = ['CompOri_',num2str(ii)];
% 
% figure(99)
thetas = 0:5:90;
lw = 2.5;

%{
subplot(6,2,1)
imagesc(Mill_Results.zThetaPref)
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Miller')

subplot(6,2,2) 
histogram(Mill_Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(Mill_Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
%legend('Neighbors', 'Random')
title('Miller')

subplot(6,2,3)
imagesc(Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Dyn')

subplot(6,2,4) 
histogram(Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
%legend('Neighbors', 'Random')
title('Dyn')

subplot(6,2,5)
imagesc(ell_Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('Ell+Dyn')

subplot(6,2,6) 
histogram(ell_Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(ell_Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
%legend('Neighbors', 'Random')
title('Ell+Dyn')

subplot(6,2,7)
imagesc(rot_Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('rot+Dyn')

subplot(6,2,8) 
histogram(rot_Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(rot_Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
%legend('Neighbors', 'Random')
title('rot+Dyn')


subplot(6,2,9)
imagesc(squi_Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('squi+Dyn')

subplot(6,2,10) 
histogram(squi_Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(squi_Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
%legend('Neighbors', 'Random')
title('squi+Dyn')

subplot(6,2,11)
imagesc(squi_rot_Results.zThetaPref);
axis square
set(gca, 'Colormap', twilight)
colorbar
title('squi+rot+Dyn')

subplot(6,2,12) 
histogram(squi_rot_Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
hold on
histogram(squi_rot_Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
hold off
axis square
%legend('Neighbors', 'Random')
title('squi+rot+Dyn')
%}

x_edge = 14-1;
y_edge = 14-1;
Nx_arbor = Inp.Nx_arbor;
zoom_inds = Inp.Nx_arbor:6*Inp.Nx_arbor;
midpt = length(zoom_inds/2);

if small_comp
    
    ang_ticks = 0:45:180
    fontsize = 14
    
    figure(100)
    
    subplot(2,1,1)
    title('With cortical dynamics')
    subplot(2,1,2)
    title('Without cortical dynamics')
    
    subplot(2,2,3)
    imagesc(Mill_Results.zThetaPref)
    axis square
    set(gca, 'Colormap', twilight)
    hold on
    rectangle('Position',[x_edge+zoom_inds(1)/Nx_arbor, y_edge+zoom_inds(1)/Nx_arbor,  zoom_inds(end)/Nx_arbor, zoom_inds(end)/Nx_arbor],...
        'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
    hold off
    colorbar('Ticks',ang_ticks, 'TickLabels', num2cell(ang_ticks))
    axis off
    ylabel('Sans cortical dynamics', 'FontSize', fontsize)
    text(-0.15, 1, 'C', 'Units','normalized', 'FontSize', fontsize)
%     title(')
    
    subplot(2,2,4)
    histogram(Mill_Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
    hold on
    histogram(Mill_Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
    hold off
    axis square
    xlim([0, max(thetas)])
    legend({'Neighbors', 'Random'}, 'Box','off')
    %title('Miller')
    ylabel('Probability')
    xlabel('\DeltaOrientation preference ( ^o )')
    text(-0.3, 1, 'D', 'Units','normalized', 'FontSize', fontsize)
    
    subplot(2,2,1)
    imagesc(Results.zThetaPref);
    axis square
    set(gca, 'Colormap', twilight)
    hold on
    rectangle('Position',[x_edge+zoom_inds(1)/Nx_arbor, y_edge+zoom_inds(1)/Nx_arbor,  zoom_inds(end)/Nx_arbor, zoom_inds(end)/Nx_arbor],...
        'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
    hold off
    colorbar('Ticks',ang_ticks, 'TickLabels', num2cell(ang_ticks))
    axis off
    ylabel('With cortical dynamics', 'FontSize', fontsize)
    text(-0.15, 1, 'A', 'Units','normalized', 'FontSize', fontsize)
    %title('Dyn')
    
    subplot(2,2,2)
    histogram(Results.DThetaNeighbors, thetas, 'Normalization', 'probability', 'FaceColor', 'k')
    hold on
    histogram(Results.DThetaRandos, thetas, 'Normalization', 'probability','DisplayStyle', 'stairs', 'EdgeColor','red', 'LineWidth',lw)
    hold off
    text(-0.3, 1, 'B', 'Units','normalized', 'FontSize', fontsize)
    axis square
    legend({'Neighbors', 'Random'}, 'Box','off')
    xlim([0,max(thetas)])
    ylabel('Probability')
    xlabel('\DeltaOrientation preference ( ^o )' )%, 'Interpreter','latex')
    %title('Dyn')
    

end