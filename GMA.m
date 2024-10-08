function [strout] = GMA(Me,Ce,Ke,tau,fls)
% Function for the execution of the Generalised Modal Analysis and for the  
% evaluation of the modal correlation coefficients
%
% Input:
% Me, Ce, Ke: (nxn) mass, damping and stiffness matrices, respectively;
% tau :       (nx1) location vector;
% Rt :        (sxn) tranformation matrix for the definition of the output quantities of interest;
% fls :        flag sismica [0] la forzante � applicata ai nodi; [1] la forzante � sismica  (capire se mantenere)
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

% common parameters
n=size(Me,2); 
% system matrices in the state-variable form 
A=[-Ke, 0*eye(n) ; 0*eye(n) , Me];
B=[0*eye(n), -Ke ; -Ke, Ce];
% location vector 
if fls==0, V=[zeros(n,1) ; tau]; end
if fls==1, V=[0*eye(n) ; eye(n)]*Me*tau; end
G=A\V;  
D=A\B;
% Generalised modal analysis
[Fi,Gam]=eig(-D);
P=Fi\G;
w0=abs(diag(Gam));
z=-real(diag(Gam))./w0;
% sorting of modal quantities
ib=1:2:2*n;  % to be automated
Gs=diag(Gam); Gs=Gs(ib);
ws=w0(ib); zs=z(ib); wds=ws.*sqrt(1-zs.^2);
Fs=Fi(1:n,ib);
Ps=P(ib);
% sorting of the modes in ascending order 
[ws,is]=sort(ws,'ascend');
zs=zs(is); Gs=Gs(is); wds=wds(is);
Ps=-Ps(is); Fs=Fs(:,is);
% building of the output data structure 
strout.Gs=Gs;
strout.Fs=Fs;
strout.ws=ws;
strout.zs=zs;
strout.wds=wds;
strout.Ps=Ps;
end

