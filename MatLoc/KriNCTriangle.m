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
% KriNCTriangle.m:
% Matrices de raideur interne EF de CR P1-NC (ordre=1) ou FS (ordre=2) 
%                             pour un triangle donne
% 
% Pour tester : aire=0.5; EdgNorm=[1,1;-1,0;0,-1]; Lga2=[2,1,1];
%
% SYNOPSIS [KNC,KCR,KFS] = KriNCTriangle(aire,EdgNorm,Lga2,ordre)
%          
% INPUT  - aire = aire du triangle
%        - EdgNorm(3,2) : coordonnees (x, y) des faces normales
%        - Lga2(3,1)    : longueurs des aretes au carre
%        - ordre            : ordre d'approximation
% OUTPUT - KNC : matrice de raideur locales pour EF de NC de l'ordre concerne
%        - KCR : matrice de raideur locales pour EF de CR P1-NC
%        - KFS : matrice de raideur locales pour EF de FS P2+BulleTriangle 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KNC,KCR,KFS]=KriNCTriangle(aire,EdgNorm,Lga2,ordre)
%
% MATRICES DE RAIDEUR
%
[KLG,K1,K2]=KriLGTriangle(aire,EdgNorm,Lga2,ordre);
KCR=4*K1;
KFS=[];
if (ordre==1)
  KNC=KCR;
end
if (ordre==2)
   Ktt=3*sum(Lga2)/(4*aire);
   K2t=zeros(6,1);
   for i=1:3        
       j=mod(i,3)+1;
       k=6-(i+j);
       K2t(i)=-0.5*KCR(i,i);
       K2t(i+3)=-KCR(j,k);
   endfor   
   KFS=[K2,K2t;K2t',Ktt];
   KNC=KFS;
   % CHECK
##   XY=[0,0;1,0;0,1];
##   [xyp,wp,lambda,np]=IntTri_Ham7(XY);
##   Glambda=-(0.5/aire)*EdgNorm;
##   wp=wp*aire;
##   GphiT=-6*lambda*Glambda;
##   ij=[2,3,1];
##   for i=1:3
##     GphiS=(4*lambda(:,i)-1)*Glambda(i,:);
##     j=ij(i); k=6-(i+j);
##     K2t(i)=wp*(sum((GphiT.*GphiS),2));
##     GphiM=4*(lambda(:,j)*Glambda(k,:)+lambda(:,k)*Glambda(j,:));
##     K2t(i+3)=wp*(sum((GphiT.*GphiM),2));
##   endfor
##   Ktt=wp*(sum((GphiT.*GphiT),2));
##   KFS2=[K2,K2t;K2t',Ktt]
endif
