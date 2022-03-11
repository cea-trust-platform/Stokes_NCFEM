%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% IntEdg_Sim4.m:
% Points d'integration sur l'arete de sommets de coordonnees XY(2,2)
% Formule d'integration sur une arete avec interpolation a 4 points 
%         exacte pour les polynomes d'ordre 3
% INPUT  - XY(2,2) : coordonnees des sommets de l'arete
% OUTPUT - xyp(4,2) : coordonnees des points d'interpolation
%        - wp(1,4) : poids associes
%        - lambda(4,2) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 4 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntEdg_Sim4(XY)

np=4;

pM=1/8;
pP=3/8;
wp=[pM,pM,pP,pP];
lambda=[0,1;
        1,0;
        1/3,1-1/3;
        1-1/3,1/3];
xyp=lambda*XY;
