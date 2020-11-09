%Caleb Holt -- May 2018
% To measure the RF spatial correlations between cells. Should include all
% neighbors out to a distance of R_arbor away. Then make a histogram of all
% the values. 

function spatialCorrRF = RFcorr(Inp, WRFon, WRFoff)

nDimV1 = Inp.nDimV1;
Nx_arbor = Inp.Nx_arbor;
R = sqrt(2)+0.1; % Inp.R_arbor/2; %radius that defines the neighborhood of cells I'm looking in. 
nV1 = Inp.nV1;
d = Inp.d;
nDimLGN = Inp.nDimLGN;
Rx_arbor = (Nx_arbor -1)/2; 

X = zeros(nDimV1^d, 1); %for periodic boundary conditions, these are bounded between 1 and L/2
Y = zeros(nDimV1^d, 1);


for i = 1:nDimV1^d
    [X(i), Y(i)] = ind2sub([nDimV1, nDimV1], i);
end


X = X-1;
Y = Y-1; 

XX = X;
YY = Y;


% if Inp.PBC
%     X(X > Inp.L) = 2*Inp.L - X(X > Inp.L); %wrap everything back into the appropriate boundaries.
%     Y(Y > Inp.L) = 2*Inp.L - Y(Y > Inp.L);
% end
% 

WRFploton = reshape(WRFon,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFploton = permute(WRFploton, [3,1,4,2]);
WRFploton = reshape(WRFploton, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

WRFplotoff = reshape(WRFoff,nDimV1,nDimV1,Nx_arbor,Nx_arbor);%all RF's
WRFplotoff = permute(WRFplotoff, [3,1,4,2]);
WRFplotoff = reshape(WRFplotoff, [nDimV1*Nx_arbor , nDimV1*Nx_arbor]);

WRF = WRFploton - WRFplotoff; % this is the actual RF looking thing. It's how Miller defines RF's in his '94 paper. 

A = reshape(Inp.A(1,:), Inp.nDimLGN, nDimLGN);
A = fftshift(A); %centers A. 
A = A((-Rx_arbor:Rx_arbor)+(round(nDimLGN/2) +1),(-Rx_arbor:Rx_arbor)+(round(nDimLGN/2) +1));%trims away cells that fall outside of the boxed arbor (i.e. Nx_arbor x Nx_arbor) 
%Find the number of pairs
% if R ~= Inp.R_arbor
%     %this will be functionalized so that I could look at only nearest
%     %neighbors for instance (R = 1), or all 8 surrounding neurons (R =
%     %sqrt(2)).
%     Af = Inp.distV1 <= R; 
%     Rprime = R; 
%     R = floor(R); %think about the case R = sqrt(2), flooring R and using it for Nx_neurons will give a 3x3 grid of cells to look at. 
%     Nx_neurons = length(-R:R); %neurons that I'm looking at right now
%     Af = fftshift(Af); %centers Af;
%     Af = Af((-R:R)+(round(nDimLGN/2) +1),(-R:R)+(round(nDimLGN/2) +1));%trims away cells that fall outside of the boxed arbor (i.e. Nx_neurons x Nx_neurons) 
% %     nDim = Inp.nDim;
% %     L = Inp.L;
% %     c_a = Inp.c_a;
% %     [A,~,~,RFinds,Nx_arbor] = Arbor2(Xdist, R,d, nDim,nDim,L, c_a);
% %     A = reshape(A(1,:), nDim, nDim);
% else 
%     Af = A((-R:R)+(round(nDimLGN/2) +1),(-R:R)+(round(nDimLGN/2) +1)); %trims away cells that fall outside of the boxed arbor (i.e. Nx_arbor x Nx_arbor) 
% end



% Af = A((-R:R)+(round(nDimLGN/2) +1),(-R:R)+(round(nDimLGN/2) +1)) ; %makes a mask for to be overlaid on each neuron. Only want to compare nonzero parts of the RF, don't include the corners where the Arbor = 0, as that should bump up the correlation I think.
% ppn = sum(Af(:)~=0); % pairs per neuron, looking over all cells that share some input from their arbor (if R = R_arbor).
% total_pairs = ppn*Inp.nV1;

%find the spatial correlation coefficient between all pairs of neurons
% K = zeros(nV1, ppn); %K = spatial correlation coefficient

spatialCorrRF = [];
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
        
        if xdist^2 + ydist^2 <= R^2 
%         if (X(i) - X(j))^2 + (Y(i) - Y(j))^2 <= R^2 %need to use the PBC X and Y here. 
            l = l + 1;
            wrf = WRF((XX(i))*Nx_arbor +(1:Nx_arbor), (YY(i))*Nx_arbor+(1:Nx_arbor)); %use the general XX and YY so that it moves around the grid appropriately.
            wrf_comp = WRF((XX(j))*Nx_arbor +(1:Nx_arbor), (YY(j))*Nx_arbor+(1:Nx_arbor));
            p_coeffs = corrcoef(wrf, wrf_comp); 
            spatialCorrRF = [spatialCorrRF; p_coeffs(2,1)]; 
            
            if p_coeffs(1,2) ~= p_coeffs(2,1)
               warning('Weird statistics thing')
            end
        end
    end
end


% for j = 1:nDimV1 %y coordinate for the cell in question
%     for i = 1:nDimV1 %x coordinate for the cell in question
%         wrf = WRF((i-1)*Nx_arbor +(1:Nx_arbor), (j-1)*Nx_arbor+(1:Nx_arbor)); %cell in question
%         ij = sub2ind([nDimV1, nDimV1], i, j); % to linearly index the first dimension of K
%         
%         for p = 1:Nx_neurons^d %Number of pairs of cells in that I'm considering now, to not count corners see if statement. number of pairs of neurons to compare the cell in question's RF to.
%             [idx, idy]= ind2sub([Nx_neurons, Nx_neurons], p); %idx = index x, idy = index y
%              
%             if Af(idx, idy) ~=0
%                 idx = idx - round(Nx_neurons/2); %to make the index periodic w.r.t. the boundary
%                 idy = idy - round(Nx_neurons/2);
%                 
%                 x = i + idx - 1; %the -1 is there cause MATLAB starts counting at 1, not zero
%                 y = j + idy - 1; 
%                 
%                 if x < 0 
%                     x = nDimV1 + x; %wraps the index around to the corresponding positive number 
%                 elseif x >= nDimV1
%                     x = x - nDimV1;
%                 end
%                 if y < 0 
%                     y = nDimV1 + y; %wraps the index around to the corresponding positive number
%                 elseif y >= nDimV1
%                     y = y - nDimV1; 
%                 end
%                 
%                 wrf_comp = WRF(x*Nx_arbor+(1:Nx_arbor), y*Nx_arbor+(1:Nx_arbor));
%                 p_coeffs = corrcoef(wrf(A~=0), wrf_comp(A~=0));
%                 
%                 K(ij, p) = p_coeffs(1,2);
%                 
%                 if p_coeffs(1,2) ~= p_coeffs(2,1)
%                     warning('Weird statistics thing') 
%                 end 
%             end
%         end
%     end
% end

% cK = K(:);
% cK = unique(cK);
% 
% figure(22)
% histogram(cK, 'Normalization', 'probability')
% title(['Correlation Coefficient between cells within R = ', num2str(Rprime)])

figure(23)
histogram(spatialCorrRF, 'Normalization', 'probability')
title(['Correlation Coefficient between cells within R = ', num2str(R)])
xlabel('Correlation Coefficient')
ylabel('Percentage')
end
    