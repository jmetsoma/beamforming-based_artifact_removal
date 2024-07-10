function Wbf=BF_type_cleaningFilt(X, L, A, alpha,gamma, Acomp, ncomp, typeReg, conds, smoothWithL)

% This function creates a beamforming-based data cleaning matrix, which can be used to filter out artifacts from EEG data provided that the topographies of the artifacts are known
%Input:
%
%X - input data (3D: channels x times x trials)
%L - lfm (2D: channels x dipoles), or [] if you do not use LFM
%A - topographies of removed artifacts
%alpha - choose in the interval of [0, 1]. For completely model-driven, use
% 1, and for completely data-driven, use 0.
%gamma - tikhonov regularization parameter, eg. 0.05, 0.1, needed in
%regularization type 1
%Acomp - complementary subspace used in regularization type 3. This is an
%orthogonal subspace to A, and it includes only neural data
%ncomp - number of components, needed in regularization types 2 and 3
%typeReg - regularization type: 1 - tikhonov, 2 - svd, 3 - svd in complementary
%subspace
%conds - [] if only one codition. If there are many conditions, then conds
%are the indices of trials which start the new condition
%smoothWithL - heuristic smoothing by compressing out smallest eigen values
%(as many values are included as smoothWithL)
%
% Output:
%
% Wbf - The cleaning filter. Multiply your to-be-cleaned data with this
% matrix from the left to complete the cleaning of the artifacts (whose
% topographies are A)
%
% .........................................................................
% 08 January 2022 : Johanna Metsomaa, NBE, Aalto University  
% .........................................................................
%

[C, T, R]=size(X);
if R>1
    if isempty(conds)
   X=X-mean(X,3)*1; 
    else
        for k=1:length(conds)
            if k<length(conds)
               X(:,:,conds(k):(conds(k+1)-1))=X(:,:,conds(k):(conds(k+1)-1))-mean(X(:,:,conds(k):(conds(k+1)-1)),3);
            else
               X(:,:,conds(k):R)=X(:,:,conds(k):R)-mean(X(:,:,conds(k):R),3);
            end
        end
    end
end
X2=reshape(X, C, []);

sampleCov=X2*X2'/(T*R);
sampleCov=sampleCov./trace(sampleCov)*C;

if size(L,1)~= size(X,1)
    lfmCov=0;
else
    lfmCov=L*L'; lfmCov=lfmCov./trace(lfmCov)*C;
end

covTot=sampleCov*(1-alpha)+lfmCov*alpha;
switch typeReg
    case 1 %tikhonov
        invCov=inv(covTot+eye(C)*gamma);
        Wbf=(A'*invCov*A)\(A'*invCov);
        
Wbf=eye(C)-A*Wbf;
    case 2 %svd truncation
        [u,d,~]=svd(covTot);
        d=diag(d);
       
        invCov=diag(1./d(1:ncomp));
        Ac=u(:,1:ncomp)'*A;
        Wbf=((Ac'*invCov*Ac)\(Ac'*invCov))*u(:,1:ncomp)';
        A=u(:,1:ncomp)*Ac;
        Wbf=u(:,1:ncomp)*u(:,1:ncomp)'-A*Wbf;
    case 3

        if isempty(Acomp)
            [Acomp,~,~]=svds(eye(size(X,1))-A*A', size(X,1)-size(A,2));
        end

        [u,d,~]=svd(Acomp'*covTot*Acomp);
        d=diag(d);
        
        if ncomp>0
            Wbf=(A'-((A'*covTot*Acomp)*(u(:,1:ncomp)*diag(1./d(1:ncomp))*u(:,1:ncomp)'))*Acomp');
        elseif ncomp==0
            Wbf=A';
            
        end
        
Wbf=eye(C)-A*Wbf;

end

if smoothWithL
    [u,~,~]=svds(covTot, smoothWithL);
    Wbf=u*u'*Wbf;
end

