function Momjk = SMM(strPSD,strout,flcr)
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
% Momjk:  (5xnxn) cross spectral moments hypermatrix, containing the entries Momjk(m,j,k)
%         that is the cross spectral moment of order (m-1) between the j-th 
%         and k-th modal oscillators 
%
% NB: in this revision the function only computes the real part of the even-order  
% spectral moments and the imaginary part of the odd-order spectral moments, 
% since these elements are useful to calculate the even-order cross spectral moments
% whose physical meaning is related to the varince of the displacements and velocieites
% in the geometrical space.

% common parameters
ws=strout.ws; wds=strout.wds; zs=strout.zs; 
N=length(ws);
if flcr==1; end
al=zeros(5,N,N); be=al; ga=al; de=al; lam_ii=zeros(2,N); Momjk=zeros(5,N,N);
Kjkmat=zeros(N,N); 
wx=strPSD.wx;
ex=strPSD.ex;
G0=strPSD.Gm;

for j=1:N
    % calculation of modal direct spectral moments
    % in this revision we calculate only m=0 e m=2 moments since they allow for 
    % the estimaton of the variances of the response in terms of displacements and velocity 
    %
    % determinaton of the branch
    if ws(j)<wx(1), kb=1; end
    if (ws(j)>=wx(1)&&ws(j)<wx(2)), kb=2; end
    if (ws(j)>=wx(2)&&ws(j)<wx(3)), kb=3; end
    if (ws(j)>=wx(3)), kb=4; end
    for im=1:2
        switch im
            case 1,m=0; 
            case 2,m=2;
        end
        lam_ii(im,j)=PSDSC(ws(j),G0,ex,wx)/(4*zs(j)*ws(j)^(3-m))*fi_mk(m,kb,zs(j),ex);
        for i=1:(kb-1)
            lam_ii(im,j)=lam_ii(im,j)+(1/ws(j)^4)*PSDSC(wx(i),G0,ex,wx)*wx(i)^(m+1)*gamma_mk(m,i,ex);
        end
    end
    % computation of the coefficient matrices for the evaluation of the per
    % cross modal spectral moments
    for k=1:N
        Kjk=(ws(j)^2-ws(k)^2)^2+4*zs(j)*zs(k)*(ws(j)^2+ws(k)^2)*ws(j)*ws(k)+4*(zs(j)^2+zs(k)^2)*ws(j)^2*ws(k)^2;
        Kjkmat(j,k)=Kjk;
        for r=1:5
            if r==1 % level m=0                
                al(1,j,k)=4*(zs(j)*ws(j)+zs(k)*ws(k))/Kjk;
                be(1,j,k)=2*(ws(j)^2-ws(k)^2+2*zs(j)*zs(k)*ws(j)*ws(k)+2*zs(k)^2*ws(k)^2)/(ws(k)*Kjk);
            end
            if r>1 % recursive calculation
                al(r,j,k)=-zs(k)*ws(k)*al(r-1,j,k)+wds(k)*be(r-1,j,k);
                be(r,j,k)=-zs(k)*ws(k)*be(r-1,j,k)-wds(k)*al(r-1,j,k);
            end
            ga(r,j,k)=al(r,j,k)*zs(k)*ws(k)+be(r,j,k)*wds(k);
            de(r,j,k)=al(r,j,k)*zs(k)*ws(k)-be(r,j,k)*wds(k);
        end
    end
end
% evaluation of the cross modal spectral moments
for j=1:N
    for k=1:N
        for r=1:5
            switch mod(r,2)
                case 1 % Re[lam_0,jk], Re[lam_2,jk], Re[lam_4,jk]
                    Momjk(r,j,k)=(-1)^((r-1)/2)*(lam_ii(1,j)*ga(r,k,j)*ws(j)^2+lam_ii(2,j)*de(r,k,j)+lam_ii(1,k)*ga(r,j,k)*ws(k)^2+lam_ii(2,k)*de(r,j,k))/2;
                case 0 % Im[lam_1,jk], Im[lam_3,jk]
                    Momjk(r,j,k)=(-1)^((r-2)/2)*(lam_ii(1,j)*ga(r,k,j)*ws(j)^2+lam_ii(2,j)*de(r,k,j)-lam_ii(1,k)*ga(r,j,k)*ws(k)^2-lam_ii(2,k)*de(r,j,k))/2;
            end
        end
    end
end


end

