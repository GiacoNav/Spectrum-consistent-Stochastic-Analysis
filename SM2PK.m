function pk = SM2PK(ly,p,ts)
% Computation of the of the p-th percentile of the peak value distribution 
% 
% Input:
% ly:     (sx5) matrix of the spectral moments of the quantity of interest 
% p:      percentile for the desired peak value estimation (p=0.5 -> peak mean value)
% ts:     duration (in s) of the pseudo-stationary time window
% 
% Output:
% pk:     (sx1) vector of the p-th percentile of the peak value for the s quantities of interest

% common parameters
s=length(ly);
pk=zeros(s,1);
% nu function
for i=1:s
    nk=1./(2.*pi)*sqrt(ly(i,3)/ly(i,1));
    qk=sqrt(1-(ly(i,2)^2)/(ly(i,1)*ly(i,3)));
    nuk= ts*nk/(-log(p));
    eta=sqrt(2.*log(2.*nuk*(1.-exp(-qk^(1.2)*sqrt( pi*log(2.*nuk))))));
    pk(i)=sqrt(ly(i,1))*eta;
end