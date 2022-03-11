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
% ErreursStokesSinusNC.m:
%
% Resout le probleme de Stokes "Sinus" 
%   avec les elements finis non conformes CR ou FS 

%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= npi*[ sin(npi*y)*{(b-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                                sin(npi*x)*{(b+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=sin(npi*x)*sin(npi*y);
%
% SYNOPSIS [Eu0,Eu1,Ep0] = ErreursStokesSinusNC(Ku,Mu,Mp,Uxh,Uyh,Ph);
%
% GLOBAL  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : numeros des sommets des triangles
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2      
% INPUT  - Ku(Nu,Nu) matrice de raideur de la vitesse, composantes x ou y
%        - Mp : matrice de masse pression
%        - Uxh  : vitesse approchee, composante x
%        - Uyh  : vitesse approchee, composante x
%        - Ph : pression approchee
%
% OUTPUT - Eu0 : Erreur L2 vitesse normalisee, calcul decompose
%        - Eu1 : Erreur H1 vitesse normalisee, calcul decompose
%        - Ep0 : Erreur L2 pression normalisee, calcul decompose
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eu0,Eu1,Ep0,Ehu0,Ehu1,Ehp0] = ErreursStokesSinusNC(Ku,Mu,Mp,Uxh,Uyh,Ph)
% 
global Nbpt CoorNeu RefNeu
global Nbtri Aires CoorBary NumTri TriEdg
global Nbedg CoorMil EdgTri RefEdg EdgNorm
global npi ordre
% Normes L2 de Uexact et grad(Uexact)
sinnpi=sin(npi)/npi; sin2npi=sin(2*npi)/(2*npi);
GUex2=2*npi^2*(0.25*(1-sin2npi)^2+...
      0.5*((1+sin2npi)^2)*(1.5-2*sinnpi+0.5*sin2npi));
Uex2=(1-sin2npi)*(1.5-2*sinnpi+0.5*sin2npi);
Pex2=(0.5-sin2npi)^2;
%
Tri=linspace(1,Nbtri,Nbtri); Tri=Tri';
%
if (ordre==1)
  Ndof=Nbedg;
  NumUloc=TriEdg;
  NumPloc=Tri;
endif
%
if (ordre==2)
  NP2=Nbpt+Nbedg;
  Ndof=NP2+Nbtri; 
  NumUloc=[NumTri,Nbpt+TriEdg,Tri+NP2];
  NumPloc=linspace(1,3*Nbtri,3*Nbtri);
  NumPloc=reshape(NumPloc,3,Nbtri)';
endif
%
% Erreur ||U_h-U_exact||_0 et ||U_h-U_exact||_1
Eu0=0; Eu1=0; Ep0=0;
Ndof_u=size(Uxh,1);
GUhGUex=0; UhUex=0; PhPex=0;
NP2=Nbedg+Nbpt;
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
PhiM=zeros(7,3);
for t=1:Nbtri
  IGLO=NumTri(t,:);  CoorNeuT=CoorNeu(IGLO,:);
  AGLO=TriEdg(t,:);  EdgNormT=EdgNorm(AGLO,:);
  UGLO=NumUloc(t,:); UX=Uxh(UGLO,1); UY=Uyh(UGLO,1);
  PGLO=NumPloc(t,:); PT=Ph(PGLO,1);
  for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
        EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
  end
  EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
  %
  aire=Aires(t);
  GxLambda=-(0.5/aire)*ones(7,1)*EdgNormT(:,1)';
  GyLambda=-(0.5/aire)*ones(7,1)*EdgNormT(:,2)';
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  sinxy=sin(npi*xyp); cosxy=cos(npi*xyp);
  cosx=cosxy(:,1); cosy=cosxy(:,2); 
  sinx=sinxy(:,1); siny=sinxy(:,2);
  awp=aire*wp';
  %
  csx=awp.*(1-cosx); csy=awp.*(cosy-1);
  Uxex=csx.*siny; Uyex=csy.*sinx;
  sinxsiny=awp.*sinx.*siny;
  GxUxex=npi*sinxsiny;  GyUxex=npi*cosy.*csx;
  GxUyex=npi*cosx.*csy; GyUyex=-GxUxex; 
  Pex=sinxsiny;
  %
  % CALCULS DE (Uh,Uex)_{0,T} (Ph,pex)_{0,T} (\grad Uh,\grad Uex)_{0,T}
  if (ordre==1)
    PhiNC=1-2*lambda;
    PhiP=ones(7,1);
    GxPhiNC=-2*GxLambda;
    GyPhiNC=-2*GyLambda;
  end
  if (ordre==2)
    lambda2=lambda.*lambda;
    PhiS=2*lambda2-lambda;
    PhiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    PhiT=2-3*sum(lambda2,2);
    PhiNC=[PhiS,PhiM,PhiT];
    PhiP=lambda;
    %
    GPhiS=4*lambda-1; 
    GxPhiS=GPhiS.*GxLambda; 
    GyPhiS=GPhiS.*GyLambda;
    %
    GxPhiM(:,ILOC)=4*(lambda(:,JLOC).*GxLambda(:,KLOC)+lambda(:,KLOC).*GxLambda(:,JLOC));
    GyPhiM(:,ILOC)=4*(lambda(:,JLOC).*GyLambda(:,KLOC)+lambda(:,KLOC).*GyLambda(:,JLOC));
    GxPhiT=-6*sum(PhiP.*GxLambda,2);
    GyPhiT=-6*sum(PhiP.*GyLambda,2);
    %
    GxPhiNC=[GxPhiS,GxPhiM,GxPhiT];
    GyPhiNC=[GyPhiS,GyPhiM,GyPhiT];
  end
  UhUex+=Uxex'*(PhiNC*UX)+Uyex'*(PhiNC*UY);
  %
  GUhGUex+=GxUxex'*(GxPhiNC*UX)+GyUxex'*(GyPhiNC*UX)+GxUyex'*(GxPhiNC*UY)+GyUyex'*(GyPhiNC*UY);
  PhPex+=Pex'*(PhiP*PT);
end
% 
Eu0=sqrt((Uxh'*Mu*Uxh+Uyh'*Mu*Uyh-2*UhUex+Uex2)/Uex2);
Eu1=sqrt((Uxh'*Ku*Uxh+Uyh'*Ku*Uyh-2*GUhGUex+GUex2)/GUex2);
Ep0=sqrt((Ph'*Mp*Ph-2*PhPex+Pex2)/Pex2);
%
