function strPSD = RS2PSD(S0,a0,TT,z,ts)
% Function for the evaluation of the parameters of a PSD function consistent
% with an assigned Response Spectrum (RS)
% 
% Input:
% S0:  the ZPA value of the RS;
% a0:  the dynamic amplification factor of the RS;
% TT:  a (3x1) vector containing the periods delimitating the branches of the RS
% z:   the value of damping ratioat which the RS is defined;
% ts:  the length of the time-window in which the earthquake can be assumed as stationary (i.e. ts=20 s)
% 
% Output: 
% strPSD: data structure that defins the RS-consistent PSD function and 
%         containing the following fields:
%         Gm: the peak value of PSD function at w=2*pi/T2
%         wx: a (3x1) vector of the circular frequencies delimitating the branches of the PSD function  
%         ex: a (4x1) vector of the exponents of the PSD model

% common parameters
kkk=[1 2]; % shape factors of the RS
ww=2*pi./TT;
l0=4.*z/(pi-4.*z);

% evaluation of PSD parameters
ex=zeros(4,1);
ex(1)=2*kkk(2)-1-elle(ww(3),z,ts);
ex(2)=2*kkk(1)-1-elle(ww(2),z,ts);
u12=(ww(3)/ww(2))^(1+ex(2))*(l0+ex(1)+1)/(1+ex(1))+(1-(ww(3)/ww(2))^(1+ex(2)))*(l0+ex(2)+1)/(1+ex(2));
ex(3)=-1-l0-u12*elle(ww(2),z,ts);
u123=(ww(2)/ww(1))^(1+ex(3))*u12+(1-(ww(2)/ww(1))^(1+ex(3)))*(l0+ex(3)+1)/(ex(3)+1);
ex(4)=-1-l0-u123*(elle(ww(1),z,ts)+2*(a0-1)/a0);
Gm=l0/(ww(2)*u12)*(a0*S0/eta0(ww(2),z,ts))^2;
% building the output data structure
strPSD.Gm=Gm;
strPSD.wx=ww;
strPSD.ex=ex;
end

% function L(w)
function yy=elle(om,zita,tis)
nm=2*(1-bii(om,zita,tis))*(log(2*nvm(om,tis)))^0.5+bii(om,zita,tis)*pi^0.5*delv(zita)^1.2d0;
dn=2*(1-bii(om,zita,tis))*(log(2*nvm(om,tis)))^0.5*log(2*nvm(om,tis)*((1-bii(om,zita,tis))));
yy=nm/dn;
end 

% function B(w)
function xx=bii(om,zita,tis)
xx=exp(-delv(zita)^1.2*(pi*log(2*nvm(om,tis)))^0.5);
end

% function for the calculation of the Vanmarcke's eta coefficient
function yy=eta0(w,z,Ts)
yy=sqrt(2*log(2*nvm(w,Ts)*(1-exp(-delv(z)^1.2*sqrt(pi*log(2*nvm(w,Ts)))))));
if w<0.37, yy=0; end		
end 

% function for the approx. calculation of the Vanmarcke's delta coefficient 
function zz=delv(z)
zz=sqrt(1-1/(1-z^2)*(1-2/pi*atan(z/sqrt(1-z^2)))^2);
end

% function for the approx. evaluation of the Vanmarcke's nu coefficient nu 
function xx=nvm(w,Ts)
p=0.5;
xx=Ts/(2*pi)*w*(-log(p))^-1;
end 



