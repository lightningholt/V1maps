
% there was bug where diagonal neighbors weren't considered.
% This fixes that.

names = fieldnames(Big_Results);

for ff =1:length(names)
    [Big_Results.(names{ff}).DThetaNeighbors, Big_Results.(names{ff}).DThetaRandos, Big_Results.(names{ff}).DPhiNeighbors, Big_Results.(names{ff}).DPhiRandos] = ...
        diffPrefTheta(Inp,Big_Results.(names{ff}).zThetaPref, Big_Results.(names{ff}).zPhiPref);
    
    Big_Results.(names{ff}).spatialCorrRF = RFcorr(Inp, Big_Results.(names{ff}).WRFon, Big_Results.(names{ff}).WRFoff);
end
