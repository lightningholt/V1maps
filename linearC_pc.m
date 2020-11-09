%function to make linear C's in just the RF

function [Cnn, Cnn_fast, e1, e2, e3] = linearC_pc(Inp)

nThal = Inp.nThal; %
R = Inp.R_arbor; %radius of the arbor
nDim = Inp.nDimV1; 
dx = 2*Inp.L/nDim;
Nrf = length(-R:dx:R); %Number/dim in RF at the max (across the middle)
Indrf = ceil(-Nrf/2):floor(Nrf/2);

d = Inp.d;
% Neff = Nrf^d;

%calculate the distance between cells within the arbor
if d == 2
    [xs0, ys0] = meshgrid(-R:dx:R);
else
    xs0 = -R:dx:R;
end
% Nrf = size(xs0,1);

%now find the distances that correspond with the cells in the arbor
%to do that initialize the arbor
A = Inp.A;
Amask = reshape(A(1,:), nDim, nDim);
Amask2 = A > 0;

if d == 2
    A_1 = fftshift(reshape(Amask2(1,:), nDim, nDim));
    A_1 = A_1(nDim/2+Indrf+1, nDim/2+Indrf+1);
else
    A_1 = fftshift(Amask2(1,:));
    A_1 = A_1(nDim/2+Indrf+1);
end 


if d == 2
    xs = xs0(A_1);
    ys = ys0(A_1);
else
    xs = xs0(A_1);
end

Neff = length(xs);

if Neff ~= sum(Amask2(1,:))
    error('Something wrong with getting the arbor centered');
end


% N = length(xs0(:));%number of grid points
if d == 2
    Xs = repmat(xs,1,Neff);
    Ys = repmat(ys,1,Neff);

    %Now find the distances within the arbor
    arborDists = sqrt( (Xs - Xs').^2 + (Ys - Ys').^2 );% Neff x Neff matrix of distances between pairs of presynaptic neurons

else
    Xs = repmat(xs,Neff,1);
    arborDists = sqrt( (Xs - Xs').^2);
end

%Cnn = makeMexHat2(Inp.distV1, Inp.sEc, Inp.aE, Inp.sIc, Inp.aI, d, 1);% last digit 1 means don't normalize (MILL94) 
Cnn = makeMexHat2(arborDists, Inp.sEc, Inp.aE, Inp.sIc, Inp.aI, d, 1);% last digit 1 means don't normalize (MILL94) 

% Cnn[arborDists ~=0] = Inp.as*Cnn[arborDists ~= 0]; 
k = 2*pi;
% Cnn = 2*pi*besselj(0, sqrt(k*arborDists));
% Cnn_fast = makeMexHat2(arborDists, Inp.sE2, Inp.aE, Inp.sI2, Inp.aI, d, 1);
%Cnn_fast = makeMexHat2(Inp.distV1, Inp.sEc2, Inp.aE, Inp.sIc2, Inp.aI, d, 1);
Cnn_fast = makeMexHat2(arborDists, Inp.sEc2, Inp.aE, Inp.sIc2, Inp.aI, d, 1);

% Cnn_fast = 2*pi*besselj(1, sqrt(0.1*k*Inp.distV1.^2));

%{
Cnn = Amask .* Cnn;
Cnn_fast = Amask .* Cnn_fast;


%fftshift will put the first neurons arbor in the middle of the grid
A_1 = fftshift(reshape(Amask2(1,:), nDim, nDim));
A_1 = A_1(nDim/2+Indrf+1, nDim/2+Indrf+1);
%double check

figure(54)
imagesc(A_1);

% zero2D = xs0(:) * 0;

Cnn2 = fftshift(Cnn);
Cnn2 = Cnn2(nDim/2+Indrf+1, nDim/2+Indrf+1);
Cnn2 = Cnn2(A_1);

Cnn_fast2 = fftshift(Cnn_fast);
Cnn_fast2 = Cnn_fast2(nDim/2+Indrf+1, nDim/2+Indrf+1);
Cnn_fast2 = Cnn_fast2(A_1);
%}

dC_slow = 1.5*Cnn; %2;
dC_fast = 1.5*Cnn_fast;%2;

[e1, v1] = svds(dC_slow, 1);
[e2, v2] = svds(dC_fast, 1);

Cnn = makeMexHat2(Inp.distV1, Inp.sEc, Inp.aE, Inp.sIc, Inp.aI, d, 1);% last digit 1 means don't normalize (MILL94) 
Cnn_fast = makeMexHat2(Inp.distV1, Inp.sEc2, Inp.aE, Inp.sIc2, Inp.aI, d, 1);

if v1 ~= 1
%     Cnn = (1/v1) * Cnn;
%     Cnn2 = (1/v1) * Cnn2;
%     dC_slow = 1.5*Cnn2;
%     [e1, v1] = svds(dC_slow,1);
%     v1
    Cnn = 1/v1 * Cnn;
end

if v2 ~= 1
%     Cnn_fast = (1/v2) * Cnn_fast;
%     Cnn_fast2 = (1/v2) * Cnn_fast2;
%     dC_fast = 1.5*Cnn_fast2;
%     [e2, v2] = svds(dC_fast,1);
%     v2
    Cnn_fast = 1/v2 * Cnn_fast;
end

e3 = randn(length(e1),1);
e3 = e3 - (e2'*e3)*e2 - (e1'*e3)*e1;
e3 = e3/(e3'*e3);

if e2'*e1 > 0.01
    warning('e1 and e2 are not orthogonal yo')
    e2'*e1
elseif e2'*e3 > 0.01
    warning('e2 and e3 could be more orthogonal')
    e2'*e3
elseif e1'*e3 > 0.01
    warning('e1 and e3 could be more orthogonal')
    e1'*e3
end

%{
xs = xs0(A_1);
ys = ys0(A_1);

Neff = length(xs);

if Neff ~= sum(Amask(1,:))
    error('Something wrong with getting the arbor centered');
end

% N = length(xs0(:));%number of grid points

Xs = repmat(xs,1,Neff);
Ys = repmat(ys,1,Neff);

%Now find the distances within the arbor
% arborDists = sqrt( (Xs - Xs').^2 + (Ys - Ys').^2 );% Neff x Neff matrix of distances between pairs of presynaptic neurons
% arborDists = sqrt( (xs - xs').^2 + (ys - ys').^2 );% Neff x Neff matrix of distances between pairs of presynaptic neurons
true_ArborDist = sqrt( (xs0 - xs0').^2 + (ys0 - ys0').^2);
%find the correlations for those cells 
% Cnn = makeMexHat2(arborDists, Inp.sE1, Inp.aE, Inp.sI1, Inp.aI, d, 1);% last digit 1 means don't normalize (MILL94) 
Cnn = makeMexHat2(true_ArborDist, Inp.sE1, Inp.aE, Inp.sI1, Inp.aI, d, 1);% last digit 1 means don't normalize (MILL94) 
% Cnn[arborDists ~=0] = Inp.as*Cnn[arborDists ~= 0]; 
k = 2*pi;
% Cnn = 2*pi*besselj(0, sqrt(k*arborDists));
% Cnn_fast = makeMexHat2(arborDists, Inp.sE2, Inp.aE, Inp.sI2, Inp.aI, d, 1);
Cnn_fast = makeMexHat2(true_ArborDist, Inp.sE2, Inp.aE, Inp.sI2, Inp.aI, d, 1);

%this part of the code relies on Cnf = -0.5*Cnn, and likewise for the fast.
%(=> C_diff = Cnn - Cnf = 1.5*Cnn)
A_2 = fftshift(A_1);
A_2 = toeplitz(A_2(:)* 1);

Cnn = A_2 .* Cnn;
Cnn_fast = A_2 .* Cnn_fast;

dC_slow = 1.5*Cnn;
dC_fast = 1.5*Cnn_fast;

[e1, v1] = svds(dC_slow, 1);
[e2, v2] = svds(dC_fast, 1);

if v1 ~= 1
    Cnn = (1/v1) * Cnn;
    dC_slow = 1.5*Cnn;
    [e1, v1] = svds(dC_slow,1);
    v1
end

if v2 ~= 1
    Cnn_fast = (1/v2) * Cnn_fast;
    dC_fast = 1.5*Cnn_fast;
    [e2, v2] = svds(dC_fast,1);
    v2
end


if v1/v2 > 1.2
    warning('Eigenvalue 1 is too large')
    v1/v2
elseif v1/v2 <0.8
    warning('Eigenvalue 1 too small')
    v1/v2
end

figure(55)
plot(1:length(e1), e1, 1:length(e2), e2)
legend('Slow Feature', 'fast Feature')

inds_RF = find(A_2(1,:) > 0);
fullCnn = zeros(nDim, nDim);
fullCnn_fast = zeros(nDim, nDim);

fullCnn(inds_RF) = Cnn;
fullCnn_fast(inds_RF) = Cnn_fast;

fullCnn = fftshift(fullCnn);
fullCnn_fast = fftshift(fullCnn_fast);
%}
