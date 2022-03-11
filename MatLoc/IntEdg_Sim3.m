%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% IntEdg_Sim3.m:
%
% Points d'integration sur l'arete de sommets de coordonnees XY(2,2)
% Formule d'integration sur une arete avec interpolation a 3 points 
%         exacte pour les polynomes d'ordre 3
%         Simpson's 1/3 rule
% INPUT  - XY(2,2) : coordonnees des sommets de l'arete
% OUTPUT - xyp(3,2) : coordonnees des points d'interpolation
%        - wp(1,3) : poids associes
%        - lambda(3,2) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 3 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntEdg_Sim3(XY)

np=3;

pM=1/6;
pP=1-2*pM;
wp=[pM,pM,pP];
lambda=[0,1;
        1,0;
        1/2,1/2];
xyp=lambda*XY;
