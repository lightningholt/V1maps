function [PrefTheta, PrefPhi, zThetaPref, zPhiPref, zThetaSel, zPhiSel, spatialCorrRF, DThetaNeighbors,DThetaRandos, DPhiNeighbors, DPhiRandos] = extractResults(Results)

%to extract all the fields from the Results structure

% fields = fieldnames(Results); 
PrefTheta = Results.PrefTheta;
PrefPhi = Results.PrefPhi;
zThetaSel = Results.zThetaSel;
zThetaPref = Results.zThetaPref;
zPhiSel = Results.zPhiSel;
zPhiPref = Results.zPhiPref;

spatialCorrRF = Results.spatialCorrRF;

DThetaNeighbors = Results.DThetaNeighbors;
DThetaRandos = Results.DThetaRandos;
DPhiNeighbors = Results.DPhiNeighbors;
DPhiRandos = Results.DPhiRandos;

end