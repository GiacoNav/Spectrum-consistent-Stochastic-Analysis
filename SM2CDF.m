function cdf = SM2CDF(ly,pl,pu,ts,np)
% Computation of the of the cumulative distribution functon of the peak value 
% 
% Input:
% ly:     (sx5) matrix of the spectral moments of the quantity of interest 
% px:     (sx1) uper bound of the variable for the definition of the CDF
% ts:     duration (in s) of the pseudo-stationary time window
% np:     number of points at which the CDF is defined
% 
% Output:
% cdf:     (2sxnp) matrix of the CDF functions. Each block is made of two columns. 
%          the first one contains the abscissa values, the secondo one, the 
%          cumulative probability value of the peak value distribution 

% common parameters
s=length(ly);
dp=(pu-pl)/(np-1); 
cdf=zeros(np,2*s);
for i=1:s
    nk=1./(2.*pi)*sqrt(ly(i,3)/ly(i,1));
    qk=sqrt(1-(ly(i,2)^2)/(ly(i,1)*ly(i,3)));
    b=(pl(i):dp(i):pu(i))';
    cdf(:,2*i-1)=b;
    for j=1:np
        num=1.-exp(-sqrt(pi/(2*ly(i,1)))*b(j)*qk^1.2);                 
        den=exp(b(j)^2/(2.*ly(i,1)))-1.;
        ak=2.*nk*num/den;
        P0=1.-exp(-b(j)^2/ly(i,1));
        cdf(j,2*i)=P0*exp(-ak*ts);
    end
    ib=b==0; cdf(ib,2*i)=0;
end

end