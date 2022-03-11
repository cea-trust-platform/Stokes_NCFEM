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
% MassNCTriangle.m:
% Matrices de masse EF de CR P1-NC (ordre=1) ou FS (ordre=2) 
%                      pour un triangle d'aire egale a 1
%
% SYNOPSIS [MCR,MFS] = MassNCTriangle(ordre,Mass2)
%
% INPUT  - ordre : ordre d'approximation
%        - Mass2 : matrice de masse P2 (si ordre=2)
% OUTPUT - MNC : matrice de masse locale pour EF de NC de l'ordre concerne
%        - MCR : matrice de masse locale pour EF de CR P1-NC
%        - MFS : matrice de masse locale pour EF de FS P2+BulleTriangle 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MNC,MCR,MFS]=MassNCTriangle(ordre,Mass2)
%
% MATRICES DE MASSE
%
MCR=diag((1/3)*ones(3,1));
MFS=[];
if (ordre==1)
  MNC=MCR;
end
if (ordre==2)
   Mts=(-1/30)*ones(3,1);
   Mtm=1/5*ones(3,1);
   M2t=[Mts;Mtm];
   Mtt=0.4;
   MFS=[Mass2,M2t;M2t',Mtt];
   MNC=MFS;
   % Check
##   XY=[0,0;1,0;0,1];
##   [xyp,wp,lambda,np]=IntTri_Ham7(XY);
##   phiT=2-3*sum((lambda.*lambda),2);
##   ij=[2,3,1];
##   for i=1:3
##     phis=2*lambda(:,i).*lambda(:,i)-lambda(:,i);
##     Mts(i)=wp*(phis.*phiT);
##     %
##     j=ij(i);
##     k=6-(i+j);
##     phim=4*lambda(:,j).*lambda(:,k);
##     %
##     Mtm(i)=wp*(phim.*phiT);
##   endfor
##   phiT2=phiT.*phiT;
##   Mtt=wp*phiT2;
##   M2t=[Mts;Mtm];
##   MNC=[Mass2,M2t;M2t',Mtt];
end
