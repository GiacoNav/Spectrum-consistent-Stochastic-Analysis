 %% PROGRAMMA PER L'ANALISI ALEATORIA DI UNA STRUTTURA ISOLATA ALLA BASE
 clc; clear variables; close all;
 
 % l'equazione differenziale del moto di un sistema MDOIF isolato alla base 
 % eccitato al piede da un'accelerazione sismica zg è:
 % M*xdd + C*xd + K*x = F(t) = -M*l*zg

 %% Importazione parametri del sistema
 % carico le matrici del sistema
 load('Matrici_BIS.mat');
 gdl = max(size(M));
 % vettore delle incidenze
 lx=zeros(gdl,1); lx(16)=1;
 ly=zeros(gdl,1); ly(17)=1;
 % spostamenti relativi
 load('Rel_drifts.mat');
 % parametri di smorzamento
 zs=0.02; % smorzamento struttura in elevazione
 zb=0.12; % smorzamento isolatori
 ts=20.;
 %% Definizione input sismico
 % Parametri Spettro di risposta
 S0=0.248*9.80665; %m/s^2
 a0=2.411;
 z0=0.05; % smorzamento al quale è definito lo spettro 
 %k1=1; k2=2;
 TT=[0.120; 0.359; 2.592]; %s
 % Calcolo parametri PSD spettrocompatibile
 strPSD=RS2PSD(S0,a0,TT,z0,ts);
%% Analisi aleatoria sistema strutturale (analitica)
tic
stout=GMA(M,C,K,ly,1);
% degrees of freedom
Rc=eye(gdl); 
mc1=MCC(M,C,K,stout,Rc);
[mms, lm]=SMM(strPSD,stout,mc1,1);
pkx=SM2PK(lm,0.5,ts); % mean peak value
pkb=SM2PK(lm,0.01,ts); % lower bound of peak values
pku=SM2PK(lm,0.99,ts); % upper bound of peak values
cdx=SM2CDF(lm,0.9*pkb,1.1*pku,ts,200); % cumulative distribution functions
% interstory drifts
mc2=MCC(M,C,K,stout,Rd);
[mms2, lm2]=SMM(strPSD,stout,mc2,1);
pkd=SM2PK(lm2,0.5,ts);
toc
save('CDF_PRO_90_ass.txt','cdx','-ascii')
