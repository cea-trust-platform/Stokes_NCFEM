%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% IntTri_Ham7.m:
% Points d'integration dans le triangle de sommets de coordonnees XY(3,2)
% Formule d'integration dans un triangle avec interpolation a 7 points
%         exacte pour les polynomes d'ordre 5
% Biblio : formule d'Hammer-Marlow-Stroud, Dautray-Lions, tome 6, page 781
% INPUT  - XY(3,2) : coordonnees des sommets du triangle
% OUTPUT - xyp(7,2) : coordonnees des points d'interpolation
%        - wp(1,7) : poids associes
%        - lambda(7,3) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 7 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntTri_Ham7(XY)


np=7;

c=sqrt(15);
pM=(155-c)/1200;
pP=(155+c)/1200;
wp=[9./40,pM,pM,pM,pP,pP,pP];
wp(1,1)=1-sum(wp(1,2:np));
lM=(6-c)/21; lM3=1-2*lM;
lP=(6+c)/21; lP3=1-2*lP;
lambda=[1./3,1./3,1-2/3;
        lM ,lM ,lM3;
        lM ,lM3,lM ;
        lM3,lM ,lM ;
        lP ,lP ,lP3;
        lP ,lP3,lP ;
        lP3,lP ,lP ];
xyp=lambda*XY;
