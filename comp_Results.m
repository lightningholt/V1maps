function comp_Results(Big_Results, Inp, ii, dirname)

if nargin <3 
    ii = 0
    dirname ='Test/'
end


Results = Big_Results.Res;
Mill_Results =Big_Results.Mill_Res ;
ell_Results = Big_Results.ell_Res ;
rot_Results = Big_Results.rot_Res;
squi_Results = Big_Results.squi_Res ;
squi_rot_Results = Big_Results.squi_rot_Res;

figname = 'CompareHists_'+num2str(ii);

sim_bins = linspace(-1,1, 20);

figure(2121)
histogram(Mill_Results.spatialCorrRF, sim_bins);
hold on
histogram(Results.spatialCorrRF, sim_bins);
histogram(ell_Results.spatialCorrRF, sim_bins);
histogram(rot_Results.spatialCorrRF, sim_bins);
histogram(squi_Results.spatialCorrRF, sim_bins);
histogram(squi_rot_Results.spatialCorrRF, sim_bins);
hold off
legend('Miller', 'Dyn', 'Dyn+Ell', 'Dyn+rotEll', 'Dyn+squiEll', 'Dyn+squi_rotEll')


[millWRFon, millWRFoff] = readyWRFtoPlot(Mill_Results, Inp);
[WRFon, WRFoff] = readyWRFtoPlot(Results, Inp);
[ellWRFon, ellWRFoff] = readyWRFtoPlot(ell_Results, Inp);
[rotWRFon, rotWRFoff] = readyWRFtoPlot(rot_Results,Inp);
[squiWRFon, squiWRFoff] = readyWRFtoPlot(squi_Results, Inp);
[squi_rotWRFon, squi_rotWRFoff] = readyWRFtoPlot(squi_rot_Results,Inp);


x_edge = 1-1;
y_edge = 12-1;
Nx_arbor = Inp.Nx_arbor;
zoom_inds = Nx_arbor:6*Nx_arbor;

%+=========================================================================
figure(2122)
colormap('gray')

subplot(2,3, 1)
imagesc(millWRFon(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - millWRFoff(y_edge*Nx_arbor + zoom_inds, x_edge*Nx_arbor + zoom_inds))
title('Miller')

subplot(2,3,2)
imagesc(WRFon(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - WRFoff(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds))
title('Dyn')

subplot(2,3,3)
imagesc(ellWRFon(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - ellWRFoff(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds))
title('Dyn+Ell')

subplot(2,3,4)
imagesc(rotWRFon(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - rotWRFoff(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds))
title('Dyn+Rot')

subplot(2,3,5)
imagesc(squiWRFon(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - squiWRFoff(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds))
title('Dyn+squis')

subplot(2,3,6)
imagesc(squi_rotWRFon(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds) - squi_rotWRFoff(y_edge*Nx_arbor + zoom_inds, x_edge * Nx_arbor + zoom_inds))
title('Dyn+squi+rot')

%+=========================================================================
figure(2123)
colormap(redblue)

subplot(2,3, 1)
imagesc(Mill_Results.wp(:,:,1))
colorbar
axis square
title('Miller')

subplot(2,3,2)
imagesc(Results.wp(:,:,1))
colorbar
axis square
title('Dyn')

subplot(2,3,3)
imagesc(ell_Results.wp(:,:,1))
colorbar
axis square
title('Dyn+Ell')

subplot(2,3,4)
imagesc(rot_Results.wp(:,:,1))
colorbar
axis square
title('Dyn+Rot')

subplot(2,3,5)
imagesc(squi_Results.wp(:,:,1))
colorbar
axis square
title('Dyn+squis')

subplot(2,3,6)
imagesc(squi_rot_Results.wp(:,:,1))
colorbar
axis square
title('Dyn+squi+rot')

%+=========================================================================
figure(2124)
colormap(redblue)

subplot(2,3, 1)
imagesc(Mill_Results.wp(:,:,2))
colorbar
axis square
title('Miller')

subplot(2,3,2)
imagesc(Results.wp(:,:,2))
colorbar
axis square
title('Dyn')

subplot(2,3,3)
imagesc(ell_Results.wp(:,:,2))
colorbar
axis square
title('Dyn+Ell')

subplot(2,3,4)
imagesc(rot_Results.wp(:,:,2))
colorbar
axis square
title('Dyn+Rot')

subplot(2,3,5)
imagesc(squi_Results.wp(:,:,2))
colorbar
axis square
title('Dyn+squis')

subplot(2,3,6)
imagesc(squi_rot_Results.wp(:,:,2))
colorbar
axis square
title('Dyn+squi+rot')



exportpdf(2121, figname, [8,1]);
exportpdf(2122, 'RFcomp'+num2str(ii), [8,1]);
exportpdf(2123, 'SlowSel'+num2str(ii), [8,1]);
exportpdf(2124, 'FastSel'+num2str(ii), [8,1]);

