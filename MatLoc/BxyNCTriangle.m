%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% BxyNCTriangle.m:
%
% Matrices de couplage vitesse-pression EF de CR P1-NC (kT=1) ou FS (kT=2) 
%                                       pour un triangle donne
% P1NC-P0 : -\int\div\vvec q
% Bx(i,T)=-(\pa_x\phi_f,1)_{0,T}
% By(i,T)=-(\pa_y\phi_f,1)_{0,T}
%
% P2NC-P1 : \int\vvec\cdot\grad q
% Bx(i,j)=(\phi_i,\pa_x\lambda_j)_{0,T}
% By(i,j)=(\phi_i,\pa_y\lambda_j)_{0,T}  
%
% Pour tester : aire=0.5; EdgNormT=[1,1;-1,0;0,-1];
%
% SYNOPSIS [BNCx,BNCy,BCRx,BCRy,BFSx,BFSy] = BxyNCTriangle(aire,EdgNormT,kT)
%          
% INPUT  - aire          : aire du triangle
%        - EdgNormT(3,2) : coordonnees (Nx, Ny) des vecteurs "face-normale" opposes aux sommets locaux
%        - kT            : kT d'approximation
% OUTPUT - BNCx,BNCy     : matrices de couplage vitesse-pression locales pour EF de l'ordre kT concerné
%        - BCRx,BCRy     : matrices de couplage vitesse-pression locales pour EF de CR P1-NC
%        - BFSx,BFSy     : matrices de couplage vitesse-pression locales pour EF de FS P2+BulleTriangle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BNCx,BNCy,BCRx,BCRy,BFSx,BFSy]=BxyNCTriangle(aire,EdgNormT,kT)
%
% MATRICES DE COUPLAGE VITESSE-PRESSION B=[Bx,By]
%
% P1NC-P0
%
[B1x,B1y,B2x,B2y]=BxyLGTriangle(aire,EdgNormT,kT);
BCRx=-2*B1x; BNCx=BCRx;
BCRy=-2*B1y; BNCy=BCRy;
%
% MATRICES DE COUPLAGE VITESSE-PRESSION B=[Bx,By]
%
% P2NC-P1-disc
%
BFSx=[]; BFSy=[];
if (kT==2)
   BFSx=[B2x,-0.25*EdgNormT(:,1)]; BNCx=BFSx;
   BFSy=[B2y,-0.25*EdgNormT(:,2)]; BNCy=BFSy;
end