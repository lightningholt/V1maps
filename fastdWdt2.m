function dWdt = fastdWdt2(W,A,GCfft,varargin)
% function dWdt = fastdWdt(W,A,GCfft,varargin)
% Calculates the Hebbian (unnormalized and unconstrained) derivative 
% of thalamocortical weight-matrix W. 
%
% In: 
% W: thalamocortical weight-matrix of size nV1 x nLGN
% A: arbor function (a matrix the same size as W)
% GCfft: nV1 x nLGN vector: this is a sum over a bunch of (in our current case two) matrices, 
% each of which is the outer product of the fft of effective cortical interaction (linear response) 
%                          with the fft  of covariance modes at a given time-scale
% So the whole sum (i.e. GCfft) contains the contribution of all time scales.
%                          

if nargin==3
    Wf = fft2(W);
    dWdt = A.* ifft2( GCfft.* Wf );
elseif nargin==5
    nDimV1 = varargin{1};
    nDimLGN = varargin{2};
    W = reshape(W,[nDimV1,nDimV1,nDimLGN,nDimLGN]);

    Wf =  fftn(W);
    dWdt = ifftn( GCfft.* Wf );
    
    dWdt = A.* reshape(dWdt,[nDimV1^2,nDimLGN^2]);
end
