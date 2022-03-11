/****************************************************************************
* Copyright (c) 2022, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
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
