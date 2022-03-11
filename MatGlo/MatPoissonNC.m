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
% MatPoissonNC.m:
% Assemblages des matrices de masse et de raideur 
%
% avec les elements finis de Galerkin non conformes (EF de Crouzeix-Raviart ou de Fortin-Soulie)
%
% SYNOPSIS [Ku,Mu] = MatPoissonNC(Nbpt,CoorNeu,Nbtri,NumTri,NbEdg,TriEdg,EdgTri,Lga,EdgNorm,ordre,nuT)
%          
% INPUT  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(NbEdg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - Lga2(NbEdg,1)    : longueurs des aretes au carre
%        - EdgNorm(NbEdg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%        - nuT : coefficient de raideur
%
% OUTPUT - Ku : matrice de raideur
%        - Mu : matrice de masse 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu] = MatPoissonNC(CoorNeu,NumTri,TriEdg,EdgTri,Lga2,EdgNorm,Aires,ordre,nuT)
%
Nbpt =size(CoorNeu,1);
Nbtri=size(NumTri,1);
NbEdg=size(EdgTri,1);
%
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle;
[MNC,MCR,MFS]=MassNCTriangle(ordre,MT2);
if (ordre==1)
  Ndof=NbEdg;
  Numloc=TriEdg;
endif
%
if (ordre==2)
  NP2=Nbpt+NbEdg;
  Ndof=NP2+Nbtri;
  Tri=linspace(1,Nbtri,Nbtri); Tri=Tri';
  Numloc=[NumTri,Nbpt+TriEdg,Tri+NP2];
endif
%
Ku=sparse(Ndof,Ndof); Mu=sparse(Ndof,Ndof);
%
for t=1:Nbtri
  aire=Aires(t);
  AGLO=TriEdg(t,:);
  EdgNormT=EdgNorm(AGLO,:);
  Lga2T   =Lga2(AGLO,:);
  for iloc=1:2
    if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
    end
  end
  EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
  %Lga2T(3)=dot(EdgNormT(3,:),EdgNormT(3,:));
  % On suppose nu constant par element
  [KNC,KCR,KFS]=KriNCTriangle(aire,EdgNormT,Lga2T,ordre);
  KuT=nuT(t)*KNC; MuT=aire*MNC;
  %% Ku et Mu
  IGLO=Numloc(t,:);
  Ku(IGLO,IGLO)=Ku(IGLO,IGLO)+KuT;
  Mu(IGLO,IGLO)=Mu(IGLO,IGLO)+MuT;
end
