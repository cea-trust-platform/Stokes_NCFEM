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
% RHSStokesNC.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "Sinus" 
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
%
% SYNOPSIS [RHSx,RHSy]=RHSStokesSinusNC()
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles (3 numeros de sommets)
%        - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des triangles
%        - Nbedg : nombre d'aretes
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - RefEdg(Nbedg,1) : Reference de chaque arete
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - RHSx : second membre composante x
%        - RHSy : second membre composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHSx,RHSy]=RHSStokesSinusNC()
global Nbpt CoorNeu RefNeu
global Nbtri Aires NumTri CoorBary
global Nbedg TriEdg RefEdg
global nu ordre
global npi npi2
%
Nbpt=size(CoorNeu,1);
Nbtri=size(NumTri,1);
%
npi2nu=npi2*nu;
cx=npi-2*npi2nu;
cy=npi+2*npi2nu;
%
Ndof=Nbedg;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  Ndof=NP2+Nbtri;
end
RHSx=zeros(Ndof,1); RHSy=zeros(Ndof,1);
%
ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
for t=1:Nbtri
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  AGLO=TriEdg(t,:);
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  sinxy=sin(npi*xyp); cosxy=cos(npi*xyp);
  FxLOC = awp.*sinxy(:,2).*(cx*cosxy(:,1)+npi2nu); 
  FyLOC = awp.*sinxy(:,1).*(cy*cosxy(:,2)-npi2nu);
  %
  % EF de CR P1NC-P0 (ordre=1) 
  if (ordre==1)
    phiCR=1-2*lambda;
    %
    RHSx(AGLO)+=sum(FxLOC.*phiCR,1)';
    RHSy(AGLO)+=sum(FyLOC.*phiCR,1)';
    %
  end
  % EF de FS P2NC-P1 (ordre=2)
  if (ordre==2)
    lambda2=lambda.*lambda;
    phiS=2*lambda2-lambda;
    phiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    phiT=2-3*sum(lambda2,2);
    %
    RHSx(IGLO)+=sum(FxLOC.*phiS,1)';
    RHSy(IGLO)+=sum(FyLOC.*phiS,1)';
    %
    AGLO+=Nbpt;
    RHSx(AGLO)+=sum(FxLOC.*phiM,1)';
    RHSy(AGLO)+=sum(FyLOC.*phiM,1)';
    %
    TGLO=NP2+t;
    RHSx(TGLO)=sum(FxLOC.*phiT,1);  
    RHSy(TGLO)=sum(FyLOC.*phiT,1);
  end
end
%
%% Traitement des CL de Dirichlet
Dof_DiriA=find(RefEdg==0);
Dof_BordA=find(RefEdg~=0);
if (ordre==1)
  Dof_Diri=Dof_DiriA;
  Dof_Bord=Dof_BordA;
end
if (ordre==2)
  Dof_DiriS=find(RefNeu==0);
  Dof_BordS=find(RefNeu~=0);
  Dof_DiriT=linspace(1,Nbtri,Nbtri)';
  NP2=Nbpt+Nbedg;
  Dof_Diri=[Dof_DiriS;Nbpt+Dof_DiriA;...
            NP2+Dof_DiriT];
  Dof_Bord=[Dof_BordS;Nbpt+Dof_BordA];
end
Nbord=size(Dof_Bord,1);
RHSx(Dof_Bord,1)=zeros(Nbord,1);
RHSy(Dof_Bord,1)=zeros(Nbord,1);
%
