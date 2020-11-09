function [WRFploton, WRFplotoff] = readyWRFtoPlot(Results, Inp)

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
end
