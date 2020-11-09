% Arbor function

function [A,XYV1,XYLGN,RFinds,Nx_arbor] = Arbor2(dist, R, varargin)
% function [A,XYV1,XYLGN,RFinds,Nx_arbor] = Arbor(dist, R,varargin)
% dist comes from the neural distance function
% R is the radius of my arbor function
% varargin{1} = dimension (1 or 2) of retinotopic grid
% nDimV1 = varargin{2};  integer edge-size of V1 grid
% nDimLGN = varargin{3}; integer edge-size of LGN grid
% L = varargin{4}; half edge-length of retinotopic square

% Out:
% A : 0/1 arbor matrix with same size as W (W is V1 x LGN)
% XYV1 : size = [nDimV1,2] : x and y locations of all V1 grid-points in linear/vectorized order
% XYLGN : size = [2,nDimLGN] : x and y locations of all LGN grid-points in linear/row-vectorized order
% RFinds : the grand-vectorized indices of all RF's within W' (transposed W which is LGN x V1, unlike W which is V1 x LGN): 
% Convention: RF indices of the same cell are contiguous, and within this group, indices for each y=const slice of RF are contiguous


if nargin>3
    d = varargin{1};
else
    d=1;
end
   
if d==1
    A = dist <= R;
    
    %find the vectorized indices of all RF's: 
    % Convention: RF indices of the same cell are contiguous, 
    % the following is written assuming the same dx (grid-resolution) for LGN and V1 and nDimV1==nDimLGN
    if nargin>=7
        XYV1 = [];
        XYLGN = [];
        
        nDimV1 = varargin{2};
        nDimLGN = varargin{3};
        L = varargin{4};

        nV1 = nDimV1;
        nLGN = nDimLGN;
        
        dxLGN = 2*L/(nDimLGN);
        R = floor(R/dxLGN);  %actual arbor R in pixels -- what are pixels in this case? Neurons
        arbvec = (-R:R)';
        Nx_arbor = length(arbvec);%edge-length of square-RF
        RFinds = zeros(nV1 * Nx_arbor,1);
        for iV1=1:nV1 %go over all V1 locations with vectorized index
            %find LGN x and y indices of all square-RF locations:
            ixs = 1 + mod(-1 + iV1 + arbvec,nDimLGN);%bring back protruding (x,y) indices to 1:nDimLGN range, by circular shift
            RFinds((iV1-1)*Nx_arbor + (1:Nx_arbor)) = (iV1-1)*nLGN + ixs;%calculate grand-vectorized index for W' matrix (note transpose) which is LGN x V1, not V1 x LGN
        end
    end    
    
elseif d==2
    nDimV1 = varargin{2};
    nDimLGN = varargin{3};
    L = varargin{4};
    c_a = varargin{5}; %scale factor for the dendrite arbor
    
    if nargin > 7
        ANTOLIK = varargin{6};
    else
        ANTOLIK = 0;
    end
    
    rd = c_a*R;
    
    
    dxV1 = 2*L/(nDimV1);
    XV1 = -L + (0:nDimV1-1)*dxV1; %now units are mm
	YV1 = repmat(XV1,nDimV1,1);%repeated rows (so second index is y-index), note that nDimV1 appears as second argument
    XV1 = repmat(XV1',1,nDimV1);%repeated columns (so first index is x-index)
    XV1 = XV1(:);%length = nDimV1^2, different x's for the same y are consecutive
    YV1 = YV1(:);%length = nDimV1^2,
    
    dxLGN = 2*L/(nDimLGN);
    XLGN = -L + (0:nDimLGN-1)*dxLGN; %now units are mm
	YLGN = repmat(XLGN,nDimLGN,1);%repeated rows (so second index is y-index)
    XLGN = repmat(XLGN',1,nDimLGN);%repeated columns (so first index is x-index)
    XLGN = XLGN(:);%length = nDimLGN^2, different x's for the same y are consecutive
    YLGN = YLGN(:);%length = nDimLGN^2
    
    %make the LGN and V1 (x,y)-coordinates of all (unconstrained-by-arbor) synapses 
    nV1 = nDimV1^d;
    nLGN = nDimLGN^d;    
%     XV1 = repmat(XV1',nLGN,1); %size = [nLGN,nV1] with repeated rows
%     YV1 = repmat(YV1',nLGN,1); %size = [nLGN,nV1] with repeated rows
%     XLGN = repmat(XLGN,1,nV1); %size = [nLGN,nV1] with repeated columns (so columns corresponds to RF of one V1 cell)
%     YLGN = repmat(YLGN,1,nV1); %size = [nLGN,nV1] with repeated columns (so columns corresponds to RF of one V1 cell)
    XV1 = repmat(XV1,1,nLGN); %size = [nV1,nLGN] with repeated columns
    YV1 = repmat(YV1,1,nLGN); %size = [nV1,nLGN] with repeated columns
    XLGN = repmat(XLGN',nV1,1); %size = [nV1,nLGN] with repeated rows (so rows corresponds to RF of one V1 cell)
    YLGN = repmat(YLGN',nV1,1); %size = [nV1,nLGN] with repeated rows (so rows corresponds to RF of one V1 cell)
    
    Xdist = abs(XLGN-XV1);
    Xdist(Xdist > L) = 2*L - Xdist(Xdist > L);%makes X-distance periodic
    Ydist = abs(YLGN-YV1);
    Ydist(Ydist > L) = 2*L - Ydist(Ydist > L);%makes Y-distance periodic
    dist = sqrt(Xdist.^2 + Ydist.^2);
    
    if ANTOLIK
        A = (Xdist.^2 + Ydist.^2 <= R^2);%arbor mask-matrix for pill-box arbor function
    else
        A = circ_overlap(R, c_a, dist, dxV1); %Arbor function more similar to Miller94;
    end
    
%     mask = rand(size(A));
%     mask(mask < sparsity) = 1;
%     mask(mask ~= 1) = 0; 
%     
%     A = A.*mask;
    
    XYV1 = [XV1(:,1),YV1(:,1)];%different x's for the same y are consecutive
    XYLGN = [XLGN(1,:);YLGN(1,:)];%different x's for the same y are consecutive  
    
    
    %find the grand-vectorized indices of all RF's: 
    % Convention: RF indices of the same cell are contiguous, and within this group, indices for each y=const slice of RF are contiguous
    %the following is written assuming the same dx (grid-resolution) for LGN and V1 and nDimV1==nDimLGN
    R = floor(R/dxLGN);  %actual arbor R in pixels
    arbvec = (-R:R)';
    Nx_arbor = length(arbvec);%edge-length of square-RF
    %find relative x and y indices of all square-RF locations, in linearized order:
    dixs = kron(ones(Nx_arbor,1),arbvec);%length = Nx_arbor^2; all x-indices for same y are consecutive
    diys = kron(arbvec,ones(Nx_arbor,1));%length = Nx_arbor^2; all x-indices for same y are consecutive
    RFinds = zeros(nV1 * Nx_arbor^2,1);
    
    for i=1:nV1 %go over all V1 locations with vectorized index
        [ixV1 ,iyV1] = ind2sub([nDimV1,nDimV1],i);%find (x,y) integer-indices of current V1 location
        %find LGN x and y indices of all square-RF locations:
        ixys = 1 + mod(-1 + [ixV1 + dixs,iyV1 + diys],nDimLGN);%bring back protruding (x,y) indices to 1:nDimLGN range, by circular shift
        ixys_vect = sub2ind([nDimLGN,nDimLGN],ixys(:,1),ixys(:,2));%vectorize index within LGN grid
        RFinds((i-1)*Nx_arbor^2 + (1:Nx_arbor^2)) = (i-1)*nLGN + ixys_vect;%calculate grand-vectorized index for W' matrix (note transpose) which is LGN x V1, not V1 x LGN
    end
end

