function [data_sound_filt, W]=FilterOutGivenPatternsSound_BF(A, data_sound, Msound)

% This function uses beamforming to create a spatial cleaning filter which takes
% into account the smoothing step by SOUND to remove an artifact topography
% (or multiple) estimated prior to SOUND
%
% Input:
% A - artifact topographies (channels x number of artifacts)
% data_sound - data that were already cleaned by SOUND (channels, times, trials)
% Msound - SOUND cleaning matrix
%
% output:
% data_sound_filt - data after filtering out the artifacts represented by A
% W - the beamforming decomposing filter to retrieve the artifact waveform(s)
%
% .........................................................................
% 26 May 2021 : Johanna Metsomaa, BNP, University of Tübingen  
% .........................................................................
%
[C, T, R]=size(data_sound);

cov=reshape(data_sound, C, []);
cov=cov*cov'/(T*R);
[u, d, ~]=svd(cov);
d=diag(d);

figure,

bar(sqrt((d)))
title('Left-click to set the svd-regularization level')
[xn, ~, ~]=ginput(1);
xn=round(xn);
P=u(:,1:xn)';


icov=diag(1./d(1:xn));
pA=P*Msound*A;

W=inv(pA'*icov*pA)*pA'*icov*P;

data_sound_filt=data_sound-reshape((P'*pA)*(W*reshape(data_sound, C, [])),...
    C, T, R);
