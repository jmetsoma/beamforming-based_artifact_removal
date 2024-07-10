function [corrected_data, sigmas,dn, times, Mcorr] = SOUND_veryfast(data0, LFM, iter,lambda0,sigmas, estimate_just_noise)
%
% This function performs the SOUND algorithm for a given data based on beamforming filtering, which provides a very fast implementation
%
% input:
% data: channels x time points x trials (note: use single-channel
% reference, _not_ average reference)
% LFM: lead-field matrix (in same reference)
% iter: number of iterations (e.g. 10)
% lambda0: regularization factor, should be the noise-to-signal ratio, e.g.
% .1, 1, 10
%
%
% optional inputs:
% sigmas: initial guess for noise levels
% estimate_just_noise: false when you wish to have the corrected data too
%

%
% output:
% corrected data: channels x time points x trials
% dn: convergence monitoring measure. Should decrease to close to zero.
% sigmas: noise level estimates
% Mcorr: spatial correction matrix 
% times: how long did each iteration take
% .........................................................................
% 21.4.2022 Johanna Metsomaa, NBE, Aalto university
% .........................................................................

[n0, t0, r0]=size(data0);
data=reshape(data0, n0, []);

chanN = size(data,1);
if nargin < 5 || isempty(sigmas)    
    [~, sigmas] = simple_wiener_veryfast(data,1e-4, true);
end

if nargin < 6
    estimate_just_noise = 0;
end

LL = LFM*LFM';

dataCov=data*data'./size(data,2);
%sound iterations
for k=1:iter
    tic
    sigmas_old = sigmas;

   
    %BF-based spatial filter matrix (spatial filters for all channels at
    %once)
    CovL=LL + lambda0*trace(LL)/(chanN-1)*diag(sigmas.^2);
    W=inv(CovL);
    W=W*diag(1./diag(W)); %From BF
    sigmas=sqrt(sum((W'*dataCov).*W',2)); %applying filters and computing noise stds


    % Following and storing the convergence of the algorithm
    dn(k) = max(abs(sigmas_old - sigmas)./sigmas_old);
    times(k)=toc;
end


CovL=LL + lambda0*trace(LL)/(chanN-1)*diag(sigmas.^2);
Mcorr=LL/CovL;%SOUND cleaning filter matrix

if estimate_just_noise
   
    corrected_data = [];
else
    % Final data correction based on the final noise-covariance estimate.
    corrected_data = Mcorr*data;
    corrected_data=reshape(corrected_data, n0, t0, r0);
end


end


