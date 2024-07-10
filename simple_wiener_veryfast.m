function [y_solved, sigmas] = simple_wiener_veryfast(data,rf, sigmas_only)

[Nc, Nt, Nr]=size(data);
data=reshape(data, Nc, []);
C = data*data'/Nt;

c=sum(diag(C))/Nc;
C=C+c*rf*eye(Nc);
W=inv(C);
D=1./diag(W);
sigmas=sqrt(D);

if sigmas_only
    y_solved=[];
else  
    y_solved=(diag(D)*W')*data;
    y_solved=reshape(y_solved, [Nc, Nt, Nr]);
end

%for comparison, channel by channel evaluation
%if 0
% for i=1:Nc
%     ir=setdiff(1:Nc, i);
%     w(ir)=-C(i,ir)*inv(C(ir,ir));
%     w(i)=1;
%     sigmas2(i)=sqrt(w*C*w');
% end
% end


