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
% DirectSolver.m:
%
% Algorithme de resolution du systeme mixte
% ( \nu\,A     0   Bx^T)(Uxh)  (Fxh)
% (     0  \nu\,A  By^T)(Uyh) =(Fyh)
% (   -Bx   -By     Sp )(Ph )  ( 0 )
%
% SYNOPSIS [Uxh,Uyh,Ph]=DirectSolver(Kri,Bx,By,MassP,StabP,UnP,RHSx,RHSy)
%
% INPUT : - Kri(ndof,ndof)    : matrice de raideur interne (une composante)
%         - Bx(ndofP,ndof)    : matrice de couplage vitesse x-pression
%         - By(ndofP,ndof)    : matrice de couplage vitesse y-pression
%         - MassP(ndofP,ndofP): matrice de masse pression
%         - UnP(ndofP,1)      : representation discrete de la fonction 1 sur \Omega.
%         - RHSx(ndofP,ndof)  : second membre vitesse
%         - RHSy(ndofP,ndof)  : second membre pression
%         - DofU(:,1)    : numeros des degres de liberte non elimines pour Ux et Uy
%         - DofP(:,1)    : numeros des degres de liberte non elimines pour P
%         - Vol               : aire du domaine      
% OUTPUT: - Uxh(ndofU,1)      : solution vitesse x
%         - Uyh(ndofU,1)      : solution vitesse y
%         - Ph(ndofP,1)       : solution pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxh,Uyh,Ph]=DirectSolver(Kri0,Bx0,By0,MassP0,UnP0,RHSx0,RHSy0,DofU,DofP,Vol)
%
global nu
%
id = tic;
% On va calculer \nu\uvec
unsurnu=1/nu;
ndofU0=size(RHSx0,1);
Uxh=zeros(ndofU0,1); Uyh=Uxh;
ndofP0=size(MassP0,1);
sizeDofU=size(DofU,1);
sizeDofP=size(DofP,1);
Kri=Kri0; 
Bx=Bx0; By=By0; 
RHSx=RHSx0; RHSy=RHSy0;
MassP=MassP0; UnP=UnP0;
if (sizeDofU>0)
  Kri=Kri(DofU,DofU);
  RHSx=RHSx(DofU,1); RHSy=RHSy(DofU,1);
  Bx=Bx(:,DofU); By=By(:,DofU);   
end
if (sizeDofP>0)
  MassP=MassP(DofP,DofP); UnP=UnP(DofP,1);
  Bx=Bx(DofP,:); By=By(DofP,:);    
end
ndofP=size(Bx,1);
%%%%%%%%%%%%%%%%
% Renumerotation 
%%%%%%%%%%%%%%%%
sU = symrcm(Kri);
Kri_s=Kri(sU,sU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumerotation "espace-composante"
%             en "composante-espace"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndofU=size(RHSx,1);
ndofVit=2*ndofU;
Projx=sparse(ndofU,ndofVit);
Projy=sparse(ndofU,ndofVit);
for i=1:ndofU
  Projx(i,2*i-1)=1;
  Projy(i,2*i  )=1;
end
Kmat=Projx'*Kri_s*Projx+Projy'*Kri_s*Projy;
Bmat=Bx(:,sU)*Projx+By(:,sU)*Projy;
Fvec=Projx'*RHSx(sU)+Projy'*RHSy(sU);
Zp=zeros(ndofP,1);
ZZp=sparse(ndofP,ndofP);
%
Mat=[Kmat , Bmat';
    -Bmat , ZZp ];
RHS=[Fvec;Zp];
Sol=Mat\RHS; 
elapsed_time=toc(id);
fprintf('Temps de resolution (Ax=b) = %7.2e s\n',elapsed_time);
Uvec_s=Sol(1:ndofVit,1);
Uvec_s=unsurnu*Uvec_s;
Pvec=Sol(ndofVit+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renumerotations inverses
%%%%%%%%%%%%%%%%%%%%%%%%%%
indU=linspace(1,ndofU,ndofU)';
Cx=2*indU-1; Cy=2*indU;
Uxh(sU)=Uvec_s(Cx); 
Uyh(sU)=Uvec_s(Cy);
if (sizeDofU>0)
  Uxh0=Uxh; Uxh=zeros(ndofU0,1);
  Uyh0=Uyh; Uyh=zeros(ndofU0,1);
  for i=1:sizeDofU
    ind=DofU(i);
    Uxh(ind)=Uxh0(i);
    Uyh(ind)=Uyh0(i);
  end
end
if (sizeDofP>0)
  Ph0=Pvec; Ph=zeros(ndofP0,1);
  Ph(DofP)=Ph0;
else
  Ph=Pvec; 
end
Ph=Ph-(UnP0'*MassP0*Ph/Vol);
