
% Absolute spatial correlation across neighboring neurons, AKA
% row-correlations of W. Have to first find the rows which correspond to
% neighboring neurons then do the pearson correlation coefficient thing

function abs_spatialCorr = absoluteSpatialCorr(Results, Inp)

nDimV1 = Inp.nDimV1;
% Nx_arbor = Inp.Nx_arbor;
R = sqrt(2)+0.5; % Inp.R_arbor/2; %radius that defines the neighborhood of cells I'm looking in. Need to make sure it's in pixels
nV1 = Inp.nV1;
d = Inp.d;
nDimLGN = Inp.nDimLGN;

W = Results.Won - Results.Woff;

X = zeros(nDimV1^d, 1); %for periodic boundary conditions, these are bounded between 1 and L/2
Y = zeros(nDimV1^d, 1);


for i = 1:nDimV1^d
    [X(i), Y(i)] = ind2sub([nDimV1, nDimV1], i);
end

XX = X;
YY = Y;

X = X-1;
Y = Y-1; 

neighbors = 8; % the neuron is next to the 8 neurons surrounding it (including diagonal)
abs_spatialCorr = [];

l = 0; 

for i = 1:nV1
    for j = (i+1):nV1
        if Inp.PBC
            if abs(X(i) - X(j)) > Inp.nDim/2
                xdist = Inp.nDim - abs(X(i) - X(j));
            else
                xdist = X(i) - X(j);
            end
            
            if abs(Y(i) - Y(j)) > Inp.nDim/2
                ydist = Inp.nDim - abs(Y(i) - Y(j));
            else
                ydist = Y(i) - Y(j);
            end
        else
            xdist = X(i) - X(j);
            ydist = Y(i) - Y(j);
            
        end
        
        if xdist^2 + ydist^2 <= R^2 %need to use the PBC X and Y here.
            pp = corrcoef(W(i,:), W(j,:));
            abs_spatialCorr = [abs_spatialCorr; pp(2,1)]; 
        end
    end
end
