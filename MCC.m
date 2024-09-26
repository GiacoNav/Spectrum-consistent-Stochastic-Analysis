function cc=MCC(Me,Ce,Ke,strout,Rt)
% Function for the execution of the Generalised Modal Analysis and for the  
% evaluation of the modal correlation coefficients
%
% Input:
% Me, Ce, Ke: (nxn) mass, damping and stiffness matrices, respectively;
% tau :       (nx1) location vector;
% Rt :        (sxn) tranformation matrix for the definition of the output quantities of interest;
% fls :        flag sismica [0] la forzante è applicata ai nodi; [1] la forzante è sismica  (capire se mantenere)
% 
% Output:
% strout: data structure containing the following fields:
%             Gs:  (nx1) vector of the square roots of the eingenvalues (complex values);
%             Fs:  (nxn) reduced modal matrix in terms of only displacements;
%             ws:  (nx1) vector of natural frequencies of modal oscillators;
%             zs:  (nx1) vector of damping of modal oscillators;
%             wds: (nx1) vectot of damped frequecies of the modal oscillators;
%             Cij: (sxnxn) hypermatrix of the C_r,ik modal combination coefficient; 
%             Dij: (sxnxn) hypermatrix of the D_r,ik modal combination coefficient; 
%             Eij: (sxnxn) hypermatrix of the E_r,ik modal combination coefficient; 

% Retrieving of the modal information
Ps=strout.Ps;
Fs=strout.Fs;
s=size(Rt,1);
Gs=strout.Gs;
n=length(Ps);
% Evaluation of the modal combinaton coefficients 
BB=Rt*Fs*diag(Ps);
AA=-2*real(BB*conj(diag(Gs)));
CC=2*real(BB);
% evaluation of the coefficient matrices 
CCC=zeros(s,n,n); DDD=CCC; EEE=CCC;  % initialization
for d=1:s
    for i=1:n
        for j=1:n
            CCC(d,i,j)=AA(d,i)*AA(d,j);
            DDD(d,i,j)=AA(d,i)*CC(d,j)-AA(d,j)*CC(d,i);
            EEE(d,i,j)=CC(d,i)*CC(d,j);
        end
    end
end
% check if the structural system is a clasically damped one
[phi,~]=eig(Me,Ke);
ll=phi'*Ce*phi;
lld=ll-diag(diag(ll));
ff=max(max(lld));
if ff<1e-6
    DDD=zeros(s,n,n);
    EEE=DDD;
end
clear phi ll lld ff

cc.Cij=CCC;
cc.Dij=DDD;
cc.Eij=EEE;
end