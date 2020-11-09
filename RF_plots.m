%function to plot RFs -- there is a NxN tiling of grids of Nx_arbor x Nx_arbor
%to show the RF of each cell in their respective grids. 

function [WRFploton, WRFplotoff] = RF_plots(WRFon, WRFoff, Inp, varargin)

nDimV1 = Inp.nDimV1;
Nx_arbor = Inp.Nx_arbor;


x_edge = 14-1;
y_edge = 14-1;
Nx_arbor = Inp.Nx_arbor;
zoom_inds = Inp.Nx_arbor:6*Inp.Nx_arbor;
midpt = length(zoom_inds/2);

WRFploton = reshape(WRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFploton = permute(WRFploton, [3,1,4,2]);
WRFploton = reshape(WRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

WRFplotoff = reshape(WRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFplotoff = permute(WRFplotoff, [3,1,4,2]);
WRFplotoff = reshape(WRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

figure(16)
imagesc(WRFploton - WRFplotoff)
colormap('gray');
colorbar('Ticks', [min(min(WRFon-WRFoff)),max(max(WRFon-WRFoff))], 'TickLabels', {'Off', 'On'})
% title('Receptive fields across cortex')
hold on
rectangle('Position',[x_edge*Nx_arbor + zoom_inds(1), y_edge * Nx_arbor + zoom_inds(1), zoom_inds(end)-Nx_arbor, zoom_inds(end)-Nx_arbor],...
    'LineWidth',2,'LineStyle','-', 'EdgeColor', 'green')
rectangle('Position',[Nx_arbor, Nx_arbor+1, Nx_arbor, Nx_arbor],...
    'LineWidth',2,'LineStyle','-', 'EdgeColor', 'yellow')
hold off

box on
xticks([])
yticks([])
end