function [Xavecorr, Tavecorr]=removeShamBF(data_corr, iR, iRsham, iRcov, lambda, GaussEnv, widthCov, widthArtifact, U)
%
% This function removes the from one condition (active TMS) the signal
% produced by the sources activated by sham condition (the second
% condition). Data from both condtions are given as input as described in
% the following.
%
% Input:
% data_corr - 3D data (channels X times X trials), including both
% conditions
% iR - indices of active condition trials
% iRsham - indices of sham condition trials
% iRcov - indices of trials which are used in estimating the covariance
% matrix (all indices)
% lambda -regularization coefficient
% GaussEnv - true = use Gaussian envelope function as weights for the covariance matrix around the time of interest. false = use square envelope (equal weights)
% widthCov - the number of time steps the covariance matrix is computed, use odd number, eg, 13
% widthArtifact - the number of time steps the artifact topographies are taken from, use an odd number, smaller than widthCov, eg, 5
% U - use []
%
% Output:
% Xavecorr - data from which sham is removed by bf
% Tavecorr -time axis indices of Xavecorr
%
% .........................................................................
% 10 July 2024 : Johanna Metsomaa, BNP, University of Tübingen  
% .........................................................................


xd=data_corr(:,:, iRcov)-1*mean(data_corr(:,:,iRcov),3);

XaveBF=(mean(data_corr(:,:,iR),3));
XaveshamBF=(mean(data_corr(:,:,iRsham),3));

tCov=floor(widthCov/2);
tArt=floor(widthArtifact/2);

t=-tCov:tCov;
if ~GaussEnv
env=ones(1,lengths(t));
else
env=exp(-(t.^2)/3.^2); 
end
env=env./sum(env);
Call=[];

%xd=xd./sqrt(mean(xd.^2,1));
for it=1:size(data_corr,2)
    C=reshape(xd(:,it,:), size(xd,1), []);
    C=C*C';
    Call(:,:,it)=C./mean(diag(C));
end

Xavecorr=zeros(size(XaveBF));

it0=1+tCov+1;

for it=1:size(data_corr,2)-tCov-1-it0
    j=it0+it-1;
    if isempty(U)
    a=XaveshamBF(:,(j-tArt):(j+tArt));
    else
        a=U;
    end
    if 1
    C=sum(Call(:,:,j-tCov:j+tCov).*repmat(permute(env,[1 3 2]),size(data_corr,1), size(data_corr,1), 1),3);
    Ci=inv(C+eye(size(C,1))*lambda);%*mean(diag(C)));
    else
        Ci=eye(size(data,1));
    end
    Save(:,j)=inv(a'*Ci*a)*(a'*Ci*XaveBF(:,j));
    Xavecorr(:,j)=XaveBF(:,j)-a*Save(:,j);
    
end
Tavecorr=(tCov+1):(size(data_corr,2)-tCov);