function MH = makeMexHat2(dists,sigE,aE,sigI,aI,varargin)
%constructs a Mexican-hat vector or matrix with the same size as distsqs
%
% In: 
% distsqs: vector or matrix of squared distances from "origin"
% sigE = sigma of excitory gaussian
% aE = amplitude of excitory gaussian
% sigI = sigma of inhibitory gaussian
% aI = amplitude of inhibitory gaussian

if nargin>5
    dim = varargin{1};
    MILL94 = varargin{2};
end

G_E = exp(-dists.^2/(2*(sigE)^2));
if dim==1
    Nor = sum(G_E(1,:));
elseif dim==2
    if ~MILL94
        Nor = sum(G_E(:));
    else 
        Nor = 1;
    end
end
G_E = G_E/Nor; %normalize G_E

MH = aE*G_E;%"excitatory" gaussian

if aI~=0 %add the "inhibitory" gaussian
    G_I = exp(-dists.^2/(2*(sigI)^2));
    if dim==1
        Nor = sum(G_I(1,:));
    elseif dim==2
        if ~MILL94
            Nor = sum(G_I(:));
        else
            Nor = 1;
        end
    end
    G_I = G_I/Nor; %normalize G_I
    MH = MH - aI*G_I; %combines the two into a Mexican Hat
end