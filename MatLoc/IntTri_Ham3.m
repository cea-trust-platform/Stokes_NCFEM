%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% IntTri_Ham3.m:
% Points d'integration dans le triangle de sommets de coordonnees XY(3,2)
% Formule d'integration dans un triangle avec interpolation a 3 points
%         exacte pour les polynomes d'ordre 2 
% Biblio : formule d'Hammer-Stroud, Dautray-Lions, tome 6, page 781
% INPUT  - XY(3,2) : coordonnees des sommets du triangle
% OUTPUT - xyp(3,2) : coordonnees des points d'interpolation
%        - wp(1,3) : poids associes
%        - lambda(3,3) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 3 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntTri_Ham3(XY)

np=3;

wp=[1/3,1/3,1/3];
wp(1,1)=1-sum(wp(1,2:np));
lambda=[1/2,1/2,0;
        1/2,0,1/2;
        0,1/2,1/2];
xyp=lambda*XY;
