%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% IntTri_Ham4.m:
% Points d'integration dans le triangle de sommets de coordonnees XY(3,2)
% Formule d'integration dans un triangle avec interpolation a 4 points
%         exacte pour les polynomes d'ordre 3 
% Biblio : formule d'Hammer-Stroud, Dautray-Lions, tome 6, page 781
% INPUT  - XY(3,2) : coordonnees des sommets du triangle
% OUTPUT - xyp(4,2) : coordonnees des points d'interpolation
%        - wp(1,4) : poids associes
%        - lambda(4,3) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 4 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntTri_Ham4(XY)

np=4;

pM=-9/16;
pP=25/48;
wp=[pM,pP,pP,pP];
lM=1/5;
lP=1-2/5;
lambda=[1/3,1/3,1-2/3;
        lM ,lM ,lP;
        lM ,lP ,lM ;
        lP ,lM ,lM ];
xyp=lambda*XY;
