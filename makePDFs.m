%Caleb Holt July 2018
%function to load data and make some figs and pdfs

function makePDFs(figname,Inp, WRFon, WRFoff, Results)

% figname = 'Sparse_2Dgrid_Norm_fft_Miller9.mat';
nameParts = strsplit(figname, '_');

switch nameParts{1}
    case 'SimDyn'
        nameNum = nameParts{end};
        saveName = strcat(nameParts{1},'_', nameNum);
    case {'MyFrame','Feedback','Sparse'}
        nameNum = strsplit(nameParts{end}, 'r');
        nameNum = strsplit(nameNum{end}, '.');
        nameNum = nameNum{1};
        saveName = strcat(nameParts{1},'_', nameNum);
% saveRecFields = strcat(saveName, '_RecFields');
% saveSelect = strcat(saveName, '_AngleSel');
% saveOri = strcat(saveName, '_OriDiff');
% savePhase = strcat(saveName, '_PhaseDiff');
% saveCorrNeighbors = strcat(saveName, '_SpatialCorr');
end

switch nameParts{1}
    case 'Sparse'
        stitle = num2str(Inp.sparsity); %really % of connections 
        stitle = strcat('% of connections = ', stitle);
    case 'Feedback'
        stitle = num2str(Inp.MvA);
        stitle = strcat('MvA = ', stitle);
    case {'MyFrame', 'SimDyn'}
        OSI = Results.OSI;
        bipolarity = Results.bipolarity;
        saturation = Results.saturation;
        
        stitle = strcat('Mean OSI = ', num2str(OSI), '; Mean Bipolarity = ', num2str(bipolarity), '; Mean Saturation = ', num2str(saturation));
end

makeAllFigs = 1;
MapsMiller(Inp, WRFon, WRFoff, Results, makeAllFigs, stitle)
exportpdf(4, saveName, [8.5,1], './');


% exportpdf(3, saveRecFields, [4, 0.7]); %creates pdf of Receptive Fields
% exportpdf(111, saveSelect, [8, 1]);
% exportpdf(45, saveOri, [8, 0.5]);
% exportpdf(46, savePhase, [8, 0.5])
% exportpdf(23, saveCorrNeighbors, [4, 0.7]);
% exportpdf(104, saveName, [4 0.75])