%neural distance calculator for various dimensions (d)

function [nDistSq, Xdist] = neurDistMill(nDim, d, L,PBC, varargin)
%d is the dimension, nDim the number of neurons per dimension, L = length of dimension;
%PBC = Periodic Boundary Conditions (1 = on, 0 = off)


nV1 = nDim^d;
spacing = 2*L/(nDim);

X = 0:nDim-1;
X = -L + X*spacing; %now units are mm

if d == 1
    X = repmat(X, nDim, 1);
    if PBC
        X = bsxfun(@minus, X, X');
        Xdist = abs(X);
        Xdist(Xdist > L) = 2*L - Xdist(Xdist > L);% -spacing;
        %Xdist(Xdist < -L) = -2*L - Xdist(Xdist < -L) + spacing;
        nDistSq = Xdist.^2;
    else
        nDistSq = (bsxfun(@minus, X, X')).^2;
    end

    
elseif d==2
    %In the 2d case we just calculate squared distances from the (1,1)-index point (=(-L,-L) location) on the grid
    %and not (as in 1d case) all the squared distances between all grid-point pairs.
    Xdist = abs(X - (-L));%calculate X-distance from the first location along x-axis
    Xdist(Xdist > L) = 2*L - Xdist(Xdist > L);%row vector
    Ydist = repmat(Xdist,nDim,1);%[nDim,nDim] matrix with repeated rows
    
    Xdist = repmat(Xdist',1,nDim);%[nDim,nDim] matrix with repeated columns
    nDistSq = Xdist.^2 + Ydist.^2;%[nDim,nDim] matrix with 2d (possibly) periodic Euclidean distance from the (1,1)-index point (=(-L,-L) location) on the grid
    if nargin > 4
        useGlobalDist = varargin{1}
        if useGlobalDist
            nDistSq = globalDist(nDistSq, nDim, d); %makes a [nDim^d, nDim^d] matrix with each column the 2d periodic euclidean distance from the neuron(column index=linear index) to every other neuron
        end
    end
    
    

    
%     [X, Y] = ind2sub([nDim, nDim], 1:nV1); %find the position (x,y) of each v1 neuron
% %     round(sqrt(nV1)), round(sqrt(nV1)) 
%     %X  = (X -1 -floor(nDim/2))*(2*L/nDim);
%     % Y  = (Y - 1  - floor(nDim/2))*(2*L/nDim);
%     X = -L + (X-1)*spacing;
%     Y = -L + (Y-1)*spacing;
%     X = repmat(X, nV1,1);
%     Y = repmat(Y, nV1,1);
%     X = bsxfun(@minus, X, X');
%     Y = bsxfun(@minus, Y, Y');
%     X = abs(X);
%     Y = abs(Y);
%     
%     if PBC
%         X(X > L) = 2*L - X(X > L);
%         Y(Y > L) = 2*L - Y(Y > L);
%     end
%     
% %     X = X/spacing;
% %     Y = Y/spacing;
%     nDist = X.^2 + Y.^2; %bsxfun(@minus, X, X').^2 + bsxfun(@minus, Y,Y').^2;
% %     nDist = sqrt(nDist);
% %     nDist = (sqrt(nDist)*spacing).^2;
end