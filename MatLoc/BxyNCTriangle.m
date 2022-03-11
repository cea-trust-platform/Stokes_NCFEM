%{
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
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
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
