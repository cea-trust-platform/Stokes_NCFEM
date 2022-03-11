%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% BxyLGTriangle.m:
%
% Matrices de couplage vitesse-pression EF de Lagrange P1 et P2 pour un triangle donne
%
% P1-P0 : -\int\div\vvec q
% Bx(i,T)=-(\pa_x\lambda_i,1)_{0,T}
% By(i,T)=-(\pa_y\lambda_i,1)_{0,T}
%
% P2-P1 : \int\vvec\cdot\grad q
% Bx(i,j)=(\phi_i,\pa_x\lambda_j)_{0,T}
% By(i,j)=(\phi_i,\pa_y\lambda_j)_{0,T}  
% 
% Pour tester : aire=0.5; EdgNormT=[1,1;-1,0;0,-1];
%
% SYNOPSIS [B1x,B1y,B2x,B2y] = BxyLGTriangle(aire,EdgNormT,kT)
%          
% INPUT  - aire          : aire du triangle
%        - EdgNormT(3,2) : coordonnees (Nx, Ny) des vecteurs "face-normale" opposes aux sommets locaux
%        - kT            : ordre d'approximation
% OUTPUT - B1x,B1y,B2x,B2y : matrices de couplage vitesse-pression locales pour EF de Lagrange P1 P2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B1x,B1y,B2x,B2y]=BxyLGTriangle(aire,EdgNormT,kT)
%
% MATRICES DE COUPLAGE VITESSE-PRESSION B=[Bx,By]
%
% P1-P0
%
B1x=0.5*EdgNormT(:,1)';
B1y=0.5*EdgNormT(:,2)';
B2x=[];B2y=[];
%
% P2-P1
%
if (kT==2)
  unsur6=1/6;
  for i=1:3
    ij=mod(i,3)+1;
    ik=mod(i+1,3)+1;
    B2x(i,i  ) = EdgNormT(i,1)*unsur6;
    B2x(i,i+3) =-B2x(i,i);
    B2y(i,i  ) = EdgNormT(i,2)*unsur6;
    B2y(i,i+3) =-B2y(i,i);
  end
  B2x(2,4) = B2x(3,3)-B2x(1,1);
  B2x(3,4) = B2x(2,2)-B2x(1,1);
  B2x(1,5) = B2x(3,3)-B2x(2,2);
  B2x(3,5) =-B2x(3,4);
  B2x(1,6) =-B2x(1,5);
  B2x(2,6) =-B2x(2,4);
  %
  B2y(2,4) = B2y(3,3)-B2y(1,1);
  B2y(3,4) = B2y(2,2)-B2y(1,1);
  B2y(1,5) = B2y(3,3)-B2y(2,2);
  B2y(3,5) =-B2y(3,4);
  B2y(1,6) =-B2y(1,5);
  B2y(2,6) =-B2y(2,4);
end