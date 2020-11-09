% function to make heterogeneous K's. See the paper 
% Distributed network interactions and their emergence in developing neocortex
% by SMith et al 2018

%Within that paper they claim that orientation selectivity corresponds
%closely to spontaneous activity of V1's correlation structure. This
%correlation structure seems to arise out of heterogeneity in local
%connections

%What I want to do as a first pass is just introduce that every cell has a
%randomly oriented ellipse, where the minor axis is half the major. 

function M = elongK(Xdist, sigE, sigI, words)
% nDimV1 = # of neurons in V1/dimension
% dx = spacing between neurons in V1
% sigE = length scale of excitatory Gaussian
% sigI = length scale of inhibitory Gaussian
% Xdist = one of the outputs of neurDist, tells the xdist from all cells
% to the (0,0) located cell. Xdist' = Ydist, and vice versa.

if nargin < 4
    words = 'ellipse'
end

d = 2;
nDimV1 = length(Xdist);
dx = mean(abs(diff(Xdist(:,1))));
kap = sigI/sigE; 

hh = 0.8 * dx/1; %smith et al use dx = 1, so the factor on the end is to rescale it
% hh is some heterogeneity parameter which determines how much the
% ellipse shifts from cell to cell
mu1 = sigE * dx/1; %smith et al use <sigma1> = 1.8

switch words
    case 'ellipse'
        % make the same ellipse kernel for all neurons
        ecc = sqrt(3)/2; %this gives a minor axis that's half the length of the major
        sigma1 = mu1 * ones(nDimV1, nDimV1);
        sigma2 = sigma1 .* sqrt(1 - ecc.^2);
        phi = zeros(nDimV1, nDimV1);
    case 'rotated'
        % same size ellipse but rotated
        ecc = sqrt(3)/2; %this gives a minor axis that's half the length of the major
        sigma1 = mu1 * ones(nDimV1, nDimV1);
        sigma2 = sigma1 .* sqrt(1 - ecc.^2);
        phi = 180*rand(nDimV1,nDimV1);
    case 'squished'
        % Same oriented ellipse but squished
        ecc = 0.13 * hh*randn(nDimV1, nDimV1) + hh; % defines the eccentricity at every location in cortical space
        ecc(ecc>1) = 2 - ecc(ecc > 1);%min(ecc,1);
        std(ecc(:))
        sigma1 = 0*hh*mu1*randn(nDimV1, nDimV1) + mu1; 
        sigma2 = sigma1 .* sqrt(1 - ecc.^2);
        phi = zeros(nDimV1, nDimV1);
    case 'rotated and squished'
        % squished and rotated
        ecc = 0.13 * hh*randn(nDimV1, nDimV1) + hh; % defines the eccentricity at every location in cortical space
        ecc(ecc>1) = 2 - ecc(ecc > 1);%min(ecc,1);
        std(ecc(:))
        sigma1 = 0*hh*mu1*randn(nDimV1, nDimV1) + mu1; 
        sigma2 = sigma1 .* sqrt(1 - ecc.^2);
        phi = 180*rand(nDimV1,nDimV1);
        
end

M = zeros(nDimV1^d, nDimV1^d);

Ydist = Xdist';
% Xdist = Xdist(:);
% Ydist = Ydist(:);

% [X, Y] = meshgrid(X,Y);

% this is the right form for elliptical distances, if phi = 90 or 0. It
% looks weird and square otherwise.
ellDist = @(Phi, sx1, sx2) (((cosd(Phi)*Xdist - sind(Phi)*Ydist)/sx1).^2 + ((sind(Phi).*Xdist + cosd(Phi).*Ydist)/sx2).^2);


for cell = 1:nDimV1^d
    
    %Need circshifts in the x and y directions. 
    [cs_x, cs_y] = ind2sub([nDimV1, nDimV1], cell);
    
    M_e = 1./(2 * pi * sigma1(cell) .* sigma2(cell)) * exp(-0.5 * ellDist(phi(cell), sigma1(cell), sigma2(cell)));
    M_i = 1./(2 * pi * kap^2 * sigma1(cell) .* sigma2(cell)) * exp(-0.5 * ellDist(phi(cell), kap*sigma1(cell), kap*sigma2(cell)));
    temp = circshift(M_e - M_i, [cs_x-1, cs_y-1]); %subtract the one so that we don't start off shifting things
    M(cell, :) = temp(:); 
    
%     R = [cosd(phi(xx)) -sind(phi(xx)); sind(phi(xx)) cosd(phi(xx))];
%     SIG = diag([1/sigma1(xx).^2, 1/sigma2(xx).^2]);
%     
%     for yy = 1:nV1
%         d2 = [X(xx); Y(xx)] - [X(yy); Y(yy)];
%         MM(xx, yy) = 1./(2*pi*sigma1(xx) .* sigma2(xx)) .* (exp(-0.5 * (R*d2).' * SIG * (R*d2)) + (-1/kap^2)* exp( (-1/(2*kap^2)).* (R*d2).' * SIG * (R*d2)));
%     end
end
% 
% pp = randi(nDimV1^d, 10);
% 
% for p = 1:length(pp)
%     MMplot = reshape(M(pp(p),:), nDimV1, nDimV1);
%     
%     figure(19)
%     contour(MMplot)
%     colorbar
%     keyboard
% end

% 
% MMplot = reshape(M, nDimV1, nDimV1, nDimV1, nDimV1);
% % MMplot = permute(MMplot, [3,1,4,2]);
% 
% for ii = 1:nDimV1
%     for jj = 1:nDimV1
%         figure(19)
%         contour(MMplot(:,:, ii, jj))
%         colorbar
%         keyboard
%     end
% end


%{
if length(K) ~= nV1
    gK = globalDist(K, nDimV1, d);
else
    gK = K;
end

heteroMask = rand(size(gK)) < 0.5;
% heteroMask = 0.5*(heteroMask + heteroMask');

gK = heteroMask.*gK;

gK = (1/2)*(gK + gK.');
%}