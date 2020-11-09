% This is to find the normalized overlap between two circles with two radii
% RA and rd (Radius Arbor, radius dendrite) seperated by a distance h. 

function A = circ_overlap(RA, c, dist, dx)
rd = c*RA; %smaller radius 

hstar = RA -rd; %Critical distance when overlap is no longer 1. 
DA = 2*RA + 1*dx; %diameter of larger circle  (the dx includes the center neuron/pixel)
A = zeros(size(dist));
A(dist < hstar) = pi*rd^2; %should be the max overlap
Amask = and(dist >= hstar, dist <= DA/2);
A(Amask) = RA.^2 * acos((dist(Amask).^2 + RA.^2 -rd.^2)./(2*dist(Amask)*RA)) + ...
    rd.^2 * acos((dist(Amask).^2 - RA.^2 +rd.^2)./(2*dist(Amask)*rd)) ...
    - 0.5*sqrt(4*dist(Amask).^2 *RA.^2 - (RA.^2 +dist(Amask).^2 - rd.^2).^2);%pi*(RA - dist(and(dist >= hstar, dist <=DA/2))).^2;
A(dist > DA/2) = 0;
A = A./max(A(:));