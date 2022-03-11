%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% MassNCTriangle.m:
% Matrices de masse EF de CR P1-NC (ordre=1) ou FS (ordre=2) 
%                      pour un triangle d'aire egale a 1
%
% SYNOPSIS [MCR,MFS] = MassNCTriangle(ordre,Mass2)
%
% INPUT  - ordre : ordre d'approximation
%        - Mass2 : matrice de masse P2 (si ordre=2)
% OUTPUT - MNC : matrice de masse locale pour EF de NC de l'ordre concerne
%        - MCR : matrice de masse locale pour EF de CR P1-NC
%        - MFS : matrice de masse locale pour EF de FS P2+BulleTriangle 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MNC,MCR,MFS]=MassNCTriangle(ordre,Mass2)
%
% MATRICES DE MASSE
%
MCR=diag((1/3)*ones(3,1));
MFS=[];
if (ordre==1)
  MNC=MCR;
end
if (ordre==2)
   Mts=(-1/30)*ones(3,1);
   Mtm=1/5*ones(3,1);
   M2t=[Mts;Mtm];
   Mtt=0.4;
   MFS=[Mass2,M2t;M2t',Mtt];
   MNC=MFS;
   % Check
##   XY=[0,0;1,0;0,1];
##   [xyp,wp,lambda,np]=IntTri_Ham7(XY);
##   phiT=2-3*sum((lambda.*lambda),2);
##   ij=[2,3,1];
##   for i=1:3
##     phis=2*lambda(:,i).*lambda(:,i)-lambda(:,i);
##     Mts(i)=wp*(phis.*phiT);
##     %
##     j=ij(i);
##     k=6-(i+j);
##     phim=4*lambda(:,j).*lambda(:,k);
##     %
##     Mtm(i)=wp*(phim.*phiT);
##   endfor
##   phiT2=phiT.*phiT;
##   Mtt=wp*phiT2;
##   M2t=[Mts;Mtm];
##   MNC=[Mass2,M2t;M2t',Mtt];
end
