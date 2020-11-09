# V1maps
Created by Caleb Holt

Contains .m or MATLAB files to study the receptive field (RF) devlopment of the visual system of carnivorous mammals under an expanded Hebbian learning rule.

Previous work, e.g. Miller 1994, studied the RF development of the primary visual cortex using a Hebbian learning rule (cells that fire together, wire together)
and found conditions under which either smooth maps or salt-and-pepper (S&P) organizations develop for all the selectivities of stimulus features. 
However, such works were unable to predict the development of both smooth maps and S&P organizations in the same cortex for different features. 
Here we expand on those old works by considering the dynamics of the cortex itself, and find a "degree of freedom" essentially, whereby different organizations 
(e.g. smooth vs salt-and-pepper) develop for different features within the same cortex.

Key Script:
C10VisMaps.m - Parent script. Calls the other scripts to simulate the network connectivity development
