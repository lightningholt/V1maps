%function to extract mean selectivity, bipolarity (how switch-y the RF's
%are between on and off), and saturation (how many cells reach Wmax). 

function [OSI, bp, sat] = OSIandBipolarity(Inp, Won, Woff, Results)
%OSI = orientation selectivity index
%bp = bipolarity
%sat = saturation 

Wmax = Inp.Wmax;

A = Inp.A;
Aon = Inp.Aon;
Aoff = Inp.Aoff;

zSel = Results.zThetaSel;
OSI = mean(zSel(:));

WD = Won(A~=0) - Woff(A~=0); %defines the Receptive field for the cortex
Wsat = Won(Aon~=0)./(Aon(Aon~=0).*Wmax) - Woff(Aoff ~= 0)./(Aoff(Aoff~=0).*Wmax); %doing the weird division because A may not be pill-box (i.e. constant) within the rec field
% Wsat = (Won(A~=0) - Woff(A ~= 0))./(A(A~=0).*Wmax);

bp = 1 - abs(mean(sign(WD))); %measure of how the rec field switches signs. 0 would mean it never switched signs. 
%bounded between 0 and 1. 

sat = mean(abs(Wsat)); %measure of how many reach their maximum values. 
%bounded between 0 and 1. 

