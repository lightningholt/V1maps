% by Caleb Holt Feb 2018
%I need a function that'll take me from distance between -L, -L neuron and
% every other neuron and will make me a distance matrix, where each row is
% the corresponding distance from the row index i neuron to every other
% neuron. 
% So I need to circularly shift the matrix found in neurDistance2 when d
% ==2 to each neuron that's in the grid and then concatenate everything
% into a global matrix. 
% The first column, is the vectorized distV1
% The next nDim-1 columns of gDist, are neurons shifted down from one another by one. 
% The next nDim columns are neurons shifted over from -L,L by one, then
% shifted down

function gDist = globalDist(distV1, nDim, d)
gDist = zeros(nDim^d, nDim^d);

for kk = 1:nDim^d
    [i, j] = ind2sub([nDim, nDim], kk);
    i = i-1;
    j = j-1;
    
    LL_dist = circshift(distV1, [i,j]);
    gDist(kk,:) = LL_dist(:);
end

% nDim = 32;
% d = 2;
% L = 10;
% PBC = 1;
% [distV1sq, Xdist] = neurDistance2(nDim, d, L, PBC);
% distV1 = sqrt(distV1sq);

