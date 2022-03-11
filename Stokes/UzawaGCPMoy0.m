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
% UzawaGCPMoy0.m:
%
% Algorithme de resolution du systeme mixte 
% ( \nu\,A     0   -Bx^T)(Uxh) =(RHS_Ux)
% (     0  \nu\,A  -By^T)(Uyh) =(RHS_Uy)
% (     Bx    By     0  )(Ph ) =(RHS_p)
%
% SYNOPSIS [Uxh,Uyh,Ph]=UzawaGCPMoy0(Kri,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,RHS_p,Pex,DofU,Vol)
%          
% Uzawa : on resout B(A^{-1})B^T\,Ph=RHS_p-B(A^{-1})RHS par un GCP, Cholesky pour inverser A
%
% INPUT : - Kri(ndofU,ndofU): matrice de raideur interne (une composante)
%         - Bx(ndofP,ndofU) : matrice de couplage vitesse x-pression
%         - By(ndofP,ndofU) : matrice de couplage vitesse y-pression
%         - Mp(ndofP,ndofP) : matrice de masse pression
%         - invMp(ndofP,1)  : inverse matrice de masse pression
%         - UnP(ndofP,1)    : representation discrete de la fonction 1 sur \Omega.
%         - RHS_Ux(ndofU,1) : second membre vitesse x
%         - RHS_Uy(ndofU,1) : second membre vitesse y 
%         - RHS_p(ndofP,1)  : second membre pression
%         - DofU(:,1)    : numeros des degres de liberte non elimines pour Ux et Uy
%         - Vol          : aire totale
% OUTPUT: - Uxh(ndofU,1) : solution vitesse x
%         - Uyh(ndofU,1) : solution vitesse y
%         - Ph(ndofP,1)  : solution pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxh,Uyh,Ph]=UzawaGCPMoy0(Kri0,Bx0,By0,Mp,invMp,UnP,RHS_Ux,RHS_Uy,RHS_p,DofU,Vol)
  global nu
  id = tic; 
  ndofU0=size(RHS_Ux,1);
  Uxh=zeros(ndofU0,1); Uyh=Uxh;
  %
  sizeDofU=size(DofU,1);
  %
  Kri=Kri0(DofU,DofU);
  RHSx=RHS_Ux(DofU,1); RHSy=RHS_Uy(DofU,1);
  Bx=Bx0(:,DofU); By=By0(:,DofU);   
  %%%%%%%%%%%%%%%%
  % Renumerotation 
  %%%%%%%%%%%%%%%%
  sU = symrcm(Kri);
  Kri_s=Kri(sU,sU);
  %%%%%%%%%%%
  % Cholesky 
  %%%%%%%%%%%
  Uchol=chol(Kri_s);
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
  Rmat=Projx'*Uchol*Projx+Projy'*Uchol*Projy;
  Bmat=Bx(:,sU)*Projx+By(:,sU)*Projy;
  Fvec=Projx'*RHSx(sU)+Projy'*RHSy(sU);
  % Vitesse initiale
  Uvec_s = Rmat\(Rmat'\Fvec);
  % Second membre pour le systeme de resolution de la pression
  % 
  FPh = nu*RHS_p+Bmat*Uvec_s;
  normFPh=sqrt(FPh'*Mp*FPh);
  % Pression initiale : Ph=Ph-(UnP'*Mp*Ph/Vol);
  Ph=invMp*FPh; Ph=Ph-(UnP'*FPh/Vol);
  %
  Ap=BKBtp(Rmat,Bmat,Ph); r=FPh-Ap;
  %
  % direction de descente initiale
  z=invMp*r; 
  d=z;
  d=d-(UnP'*r/Vol); % d=d-(UnP'*Mp*z/Vol);
  % norme du résidu au carré
  rr=r'*r;
  %
  rz=r'*z;
  %
  eps=1.e-7;
  eps2=eps^2*rr
  k=0;
  while (rr>eps2)
    Ad=BKBtp(Rmat,Bmat,d);
    alPha=rz/(Ad'*d);
    Ph=Ph+alPha*d;
    r=r-alPha*Ad;
    z=invMp*r;
    rz1=r'*z;
    beta=rz1/rz;
    d=z+beta*d;
    d=d-(UnP'*r/Vol);% d=d-(UnP'*Mp*d/Vol);
    rr=r'*r;
    rz=rz1;
    k=k+1;
  end
  Ph=Ph-(UnP'*Mp*Ph/Vol);
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Resolution de la vitesse
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  Uvec_s=(1/nu)*(Uvec_s-Rmat\(Rmat'\(Bmat'*Ph)));
  elapsed_time=toc(id);
  fprintf('Temps de resolution (Ax=b) = %7.2e s\n',elapsed_time);
  fprintf('Nombre d iterations = %i\n',k);
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Renumerotations inverses
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  indU=linspace(1,ndofU,ndofU)';
  Cx=2*indU-1; Cy=2*indU;
  Uxh(sU)=Uvec_s(Cx); 
  Uyh(sU)=Uvec_s(Cy);
  Uxh0=Uxh; Uxh=zeros(ndofU0,1);
  Uyh0=Uyh; Uyh=zeros(ndofU0,1);
  for i=1:sizeDofU
    ind=DofU(i);
    Uxh(ind)=Uxh0(i);
    Uyh(ind)=Uyh0(i);
  end
  %
end
