function [Momjk,lamy] = SMM(strPSD,strout,cc,flcr)
% Evaluation of the spectral moment matrices in modal space for 
% RS-consistent seismic load.
% 
% Input:
% strPSD: data structure that defins the RS-consistent PSD function and 
%         containing the following fields:
%             Gm: the peak value of PSD function at w=2*pi/T2
%             wx: a (3x1) vector of the circular frequencies delimitating the branches of the PSD function  
%             ex: a (4x1) vector of the exponents of the PSD model
% strout: data structure containing the following fields:
%             Gs:  (nx1) vector of the square roots of the eingenvalues (complex values);
%             Fs:  (nxn) reduced modal matrix in terms of only displacements;
%             ws:  (nx1) vector of natural frequencies of modal oscillators;
%             zs:  (nx1) vector of damping of modal oscillators;
%             wds: (nx1) vectot of damped frequecies of the modal oscillators;
%             Cij: (sxnxn) hypermatrix of the C_r,ik modal combination coefficient (not used here); 
%             Dij: (sxnxn) hypermatrix of the D_r,ik modal combination coefficient (not used here); 
%             Eij: (sxnxn) hypermatrix of the E_r,ik modal combination coefficient (not used here); 
% flcr    flag for the evaluation of only direct (flcr=0) spectral moments 
%         in modal space, or for the calculation of the cross spectral moments too (flcr=1) 
% 
% Output:
% Momjk:  (5xnxnx2) cross spectral moments hypermatrix, containing the entries Momjk(m,j,k)
%         that is the cross spectral moment of order (m-1) between the j-th 
%         and k-th modal oscillators. The last dimension is for real (1) and imaginary (2) parts.
% lamy:   (sx5) matrix of the spectral moments of the quantity of interest 
%
% NB: in this revision the function only computes the real part of the even-order  
% spectral moments and the imaginary part of the odd-order spectral moments, 
% since these elements are useful to calculate the even-order cross spectral moments
% whose physical meaning is related to the varince of the displacements and velocities
% in the geometrical space. In the next version all the element calculation
% will be implemented

% common parameters
ws=strout.ws; wds=strout.wds; zs=strout.zs; 
N=length(ws);
if flcr==1; end
al=zeros(9,N,N); be=al; ga=al; de=al; ep=al;
lam_ii=zeros(4,N); Momjk=zeros(5,N,N,2);
Kjkmat=zeros(N,N); 
wx=flipud(strPSD.wx);
ex=strPSD.ex;
G0=strPSD.Gm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:N
    % calculation of modal direct spectral moments
    %
    % determinaton of the branch
    if ws(j)<wx(1), kb=1; end
    if (ws(j)>=wx(1)&&ws(j)<wx(2)), kb=2; end
    if (ws(j)>=wx(2)&&ws(j)<wx(3)), kb=3; end
    if (ws(j)>=wx(3)), kb=4; end
    % evaluation of direct modal spectral moments by analytical formulation
    for m=0:3  %(m is the momnet order)
        lam_ii(m+1,j)=PSDSC(ws(j),G0,ex,wx)/(4*zs(j)*ws(j)^(3-m))*fi_mk(m,kb,zs(j),ex);
        for i=1:(kb-1)
            lam_ii(m+1,j)=lam_ii(m+1,j)+(1/ws(j)^4)*PSDSC(wx(i),G0,ex,wx)*wx(i)^(m+1)*gamma_mk(m,i,ex);
        end
    end
end        
 % computation of the coefficient matrices for the evaluation of the per
    % cross modal spectral moments
for j=1:N
    for k=1:N
        Kjk=(ws(j)^2-ws(k)^2)^2+4*zs(j)*zs(k)*(ws(j)^2+ws(k)^2)*ws(j)*ws(k)+4*(zs(j)^2+zs(k)^2)*ws(j)^2*ws(k)^2;
        Kjkmat(j,k)=Kjk;
        for r=1:9
            if r==1 % level m=0                
                al(1,j,k)=4*(zs(j)*ws(j)+zs(k)*ws(k))/Kjk;
                be(1,j,k)=2*(ws(j)^2-ws(k)^2+2*zs(j)*zs(k)*ws(j)*ws(k)+2*zs(k)^2*ws(k)^2)/(wds(k)*Kjk);
            end
            if r>1 % recursive calculation
                al(r,j,k)=-zs(k)*ws(k)*al(r-1,j,k)+wds(k)*be(r-1,j,k);
                be(r,j,k)=-zs(k)*ws(k)*be(r-1,j,k)-wds(k)*al(r-1,j,k);
            end
            ga(r,j,k)=al(r,j,k)*zs(k)*ws(k)+be(r,j,k)*wds(k);
            de(r,j,k)=al(r,j,k)*zs(k)*ws(k)-be(r,j,k)*wds(k);
            ep(r,j,k)=2*be(r,j,k)*zs(k)*ws(k)*wds(k)+al(r,j,k)*(ws(k)^2-2*wds(k)^2);
        end
    end
end        
       

% evaluation of the cross modal spectral moments 
for j=1:N
    for k=1:N
        for r=1:7  % max 7 in order to have moment order up to 4 (e.g. r=1 -> m=0)
            
            switch mod(r,2)
                case 1 % even-order moments (r is odd, m is even)
                    % real part
                    Momjk(r,j,k,1)=0.5*(-1)^((r-1)/2)*(lam_ii(1,j)*ga(r,k,j)*ws(j)^2+lam_ii(3,j)*de(r,k,j)+lam_ii(1,k)*ga(r,j,k)*ws(k)^2+lam_ii(3,k)*de(r,j,k));
                    % imaginary part 
                    Momjk(r,j,k,2)=0.5*(-1)^((r-1)/2)*(lam_ii(2,j)*ep(r,k,j)-lam_ii(2,k)*ep(r,j,k)+lam_ii(4,j)*al(r,j,k)-lam_ii(4,k)*al(r,k,j));
                case 0 % odd-order moments (r is even, m is odd) 
                    % real part
                    Momjk(r,j,k,1)=-0.5*(-1)^((r+2)/2)*(lam_ii(2,j)*ep(r,k,j)+lam_ii(2,k)*ep(r,j,k)+lam_ii(4,j)*al(r,k,j)+lam_ii(4,k)*al(r,j,k));
                    % imaginary part 
                    Momjk(r,j,k,2)=-0.5*(-1)^((r+2)/2)*(lam_ii(1,j)*ga(r,k,j)*ws(j)^2+lam_ii(3,j)*de(r,k,j)-lam_ii(1,k)*ga(r,j,k)*ws(k)^2-lam_ii(3,k)*de(r,j,k));
            end
        end
    end
end
 





% Evaluation of the spectral moment of the quantity of interest in the geometrical space
s=size(cc.Cij,1);
CCC=cc.Cij;
DDD=cc.Dij;
EEE=cc.Eij;
lamy=zeros(s,5);
for i=1:s
    for r=1:5
        lamy(i,r)=sum(sum(squeeze(CCC(i,:,:)).*squeeze(Momjk(r,:,:,1))))-sum(sum(squeeze(DDD(i,:,:)).*squeeze(Momjk(r+1,:,:,2))))+sum(sum(squeeze(EEE(i,:,:)).*squeeze(Momjk(r+2,:,:,1))));
    end
end

end

function gmk=gamma_mk(m,k,ex)
%Funzione gamma(m,k) per momento spettrale diretto
gmk=(ex(k+1)-ex(k))/((1+m+ex(k+1))*(1+m+ex(k)));
end

function fmk = fi_mk(m,k,z,ex)
%Funzione fi(m,k) per la definizione del momento spettrale modale diretto
ee=0;
if (m==0||m==2), ee=ex(k); end
if (m==1), ee=2+2*ex(k); end
if (m==3), ee=8+3*ex(k); end
fmk=pi-(4*z)/((m+1+ex(k))*(1+m))*ee;
end


