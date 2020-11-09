%Caleb Holt 
%Make histograms showing the difference in angle between neighboring cells
%and random cells.

function [DThetaNeighbors, DThetaRandos, DPhiNeighbors, DPhiRandos] = diffPrefTheta(Inp, zThetaPref, zPhiPref)

nDimV1 = Inp.nDimV1;
Nx_arbor = Inp.Nx_arbor;
R = sqrt(2)+0.5; % Inp.R_arbor/2; %radius that defines the neighborhood of cells I'm looking in. Need to make sure it's in pixels
nV1 = Inp.nV1;
d = Inp.d;
nDimLGN = Inp.nDimLGN;

X = zeros(nDimV1^d, 1); %for periodic boundary conditions, these are bounded between 1 and L/2
Y = zeros(nDimV1^d, 1);


for i = 1:nDimV1^d
    [X(i), Y(i)] = ind2sub([nDimV1, nDimV1], i);
end

XX = X;
YY = Y;

X = X-1;
Y = Y-1; 
% if Inp.PBC
%     X(X > Inp.L) = 2*Inp.L - X(X > Inp.L); %wrap everything back into the appropriate boundaries.
%     Y(Y > Inp.L) = 2*Inp.L - Y(Y > Inp.L);
% end


neighbors = 8; % the neuron is next to the 8 neurons surrounding it (including diagonal)
DThetaNeighbors = [];
DThetaRandos = [];
DPhiNeighbors = [];
DPhiRandos = [];

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
%         if (X(i) - X(j))^2 + (Y(i) - Y(j))^2 <= R^2 %need to use the PBC X and Y here.

            l = l +1;
%             DThetaNeighbors = [DThetaNeighbors; abs(PrefTheta(XX(i), YY(i)) - PrefTheta(XX(j), YY(j)))];
            DThetaNeighbors = [DThetaNeighbors; abs(zThetaPref(XX(i), YY(i)) - zThetaPref(XX(j), YY(j)))];
            DPhiNeighbors = [DPhiNeighbors; abs(zPhiPref(XX(i), YY(i)) - zPhiPref(XX(j), YY(j)))];
            
            randX = randi(nDimV1 , 1);
            randY = randi(nDimV1, 1); 
            
            
%             DThetaRandos = [DThetaRandos; abs(PrefTheta(XX(i), YY(i)) - PrefTheta(randX, randY))];
            DThetaRandos = [DThetaRandos; abs(zThetaPref(XX(i), YY(i)) - zThetaPref(randX, randY))];
            DPhiRandos = [DPhiRandos; abs(zPhiPref(XX(i), YY(i)) - zPhiPref(randX, randY))];
        end
    end
end

% edges = 0:10:180;
edges = 0:10:90;
DThetaNeighbors(DThetaNeighbors > 90) = 180 - DThetaNeighbors(DThetaNeighbors > 90); % 180 is the same as 0 angles different so need to find the shortest difference in preferred angle
DThetaRandos(DThetaRandos > 90) = 180 - DThetaRandos(DThetaRandos > 90);

DPhiNeighbors(DPhiNeighbors > 180) = 360 - DPhiNeighbors(DPhiNeighbors > 180);
DPhiRandos(DPhiRandos > 180) = 360 - DPhiRandos(DPhiRandos > 180);


figure(45)
subplot(2,2,1)
histogram(DThetaNeighbors,edges,  'Normalization', 'probability')
title('Histogram of Preferred Ori Difference Neighbors')
xlabel('Diff of Preferred Ori')
ylabel('Probability')



subplot(2,2,3)
histogram(DThetaRandos,edges, 'Normalization', 'probability')
title('Histogram of Preferred Ori Difference Random Pairs')
xlabel('Diff of Preferred Ori')
ylabel('Probability')
        
subplot(1,2,2)
imagesc(zThetaPref);
title('Preferred Ori Map (arg z)')
% imagesc(PrefTheta)
% title('Preferred Ori Map (Max Resp)')
colormap('jet');
colorbar;   

figure(46)
subplot(2,2,1)
histogram(DPhiNeighbors,edges,  'Normalization', 'probability')
title('Histogram of Preferred Phase Difference Neighbors')
xlabel('Diff of Preferred Ori')
ylabel('Probability')



subplot(2,2,3)
histogram(DPhiRandos,edges, 'Normalization', 'probability')
title('Histogram of Preferred Phase Difference Random Pairs')
xlabel('Diff of Preferred Ori')
ylabel('Probability')
        
subplot(1,2,2)
imagesc(zPhiPref);
title('Preferred Phase Map (arg z)')
% imagesc(PrefTheta)
% title('Preferred Ori Map (Max Resp)')
colormap('jet');
colorbar;   
end