%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% IntEdg_Boo5.m:
%
% Points d'integration sur l'arete de sommets de coordonnees XY(2,2)
% Formule d'integration sur une arete avec interpolation a 5 points 
%         exacte pour les polynomes d'ordre 6
% Biblio : Boole's rule, Abramowitz and Stegun (1972, p. 886)
% INPUT  - XY(2,2) : coordonnees des sommets de l'arete
% OUTPUT - xyp(5,2) : coordonnees des points d'interpolation
%        - wp(1,5) : poids associes
%        - lambda(5,2) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 5 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntEdg_Boo5(XY)

np=5;

p0=7/90;
pQ=16/45;
pM=1-2*(p0+pQ);
wp=[p0,p0,pQ,pQ,pM];
lambda=[0,1;
        1,0;
        3/4,1/4;
        1/4,3/4;
        1/2,1/2];
xyp=lambda*XY;
