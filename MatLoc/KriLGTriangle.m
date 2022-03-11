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
% KriLGTriangle.m:
% Matrices de raideur interne EF de Lagrange P1 et P2 pour un triangle donne
% 
% Pour tester : aire=0.5; EdgNorm=[1,1;-1,0;0,-1]; Lga2=[2,1,1]
%
% SYNOPSIS [KLG,K1,K2] = KriLGTriangle(aire,EdgNorm,Lga2,ordre)
%          
% INPUT  - aire : aire du triangle
%        - EdgNorm(3,2) : coordonnees (x, y) des faces normales
%        - Lga2(3) : longueurs des aretes au carre
%        - ordre            : ordre d'approximation
% OUTPUT - KLG, K1, K2 : matrices de raideur locales pour EF de Lagrange ordre, P1 P2
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KLG,K1,K2]=KriLGTriangle(aire,EdgNorm,Lga2,ordre)
%
% MATRICES DE RAIDEUR
%
% P1
%
% K1 = Produit Scalaire des gradients des coordonnees barycentriques
%      
%
K1=zeros(3,3);
a4=1/(4*aire);
for i=1:3
  K1(i,i)=a4*Lga2(i);
end
K1(1,2)=a4*dot(EdgNorm(1,:),EdgNorm(2,:));
K1(2,1)=K1(1,2);
K1(1,3)=-K1(1,1)-K1(1,2);
K1(3,1)=K1(1,3);
K1(2,3)=-K1(2,2)-K1(2,1);
K1(3,2)=K1(2,3);

if (ordre==1)
   KLG=K1;
   K2=[];
end
if (ordre==2)
   K2=zeros(6,6);
   for i=1:3
     K2(i,i)=K1(i,i);
     for j=i+1:3
       K2(i,j)=-K1(i,j)/3; K2(j,i)=K2(i,j);
     end
     j=mod(i,3)+1;
     k=6-(i+j);       
     K2(i+3,i+3)=8*(K1(i,i)-K1(j,k))/3;
   end
   K2(1,5)=4*K1(1,3)/3; K2(5,1)=K2(1,5);
   K2(1,6)=4*K1(1,2)/3; K2(6,1)=K2(1,6);
   K2(2,4)=4*K1(2,3)/3; K2(4,2)=K2(2,4);
   K2(2,6)=K2(1,6);     K2(6,2)=K2(2,6);
   K2(3,4)=K2(2,4);     K2(4,3)=K2(3,4);
   K2(3,5)=K2(1,5);     K2(5,3)=K2(3,5);
   K2(4,5)=8*K1(1,2)/3; K2(5,4)=K2(4,5);
   K2(4,6)=2*K2(3,5);   K2(6,4)=K2(4,6);
   K2(5,6)=2*K2(3,4);   K2(6,5)=K2(5,6);
##   verif=0;
##   if (verif==1)
##      Un = ones(6,1);
##      KLGn=K2*Un; printf(' ordre 2, max Kri*Un = %7.2e\n',max(abs(KLGn)));
##   end
   KLG=K2;
end
