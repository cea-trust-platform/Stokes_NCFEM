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
% MatStokesNC.m:
%
% Assemblages des matrices de masse, de raideur et de couplage
%             EF de CR P1-NC (ordre=1) ou FS (ordre=2) 
%
% SYNOPSIS [Ku,Mu,Mp,invMp,Bx,By] = MatStokesNC
%          
% INPUT  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - Lga2(Nbedg,1)    : longueurs des aretes au carre
%        - EdgNorm(Nbedg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - Ku : matrice de raideur vitesse
%        - Mu : matrice de masse vitesse
%        - Mp : matrice de masse pression
%        - invMp : inverse de la matrice de masse pression
%        - Bx : matrice de couplage vitesse-pression
%        - By : matrice de couplage vitesse-pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu,Mp,invMp,Bx,By] = MatStokesNC()
%
global Nbpt CoorNeu
global Nbtri Aires NumTri TriEdg
global Nbedg Lga2 EdgNorm EdgTri
global ordre
%
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle();
[MNC,MCR,MFS]=MassNCTriangle(ordre,MT2);
if (ordre==1)
  Ndof=Nbedg;
  NumUloc=TriEdg;
  Nploc=1;
  Mploc=1; invMploc=1;
endif
%
if (ordre==2)
  NP2=Nbpt+Nbedg;
  Ndof=NP2+Nbtri;
  Tri=linspace(1,Nbtri,Nbtri)';
  NumUloc=[NumTri,Nbpt+TriEdg,Tri+NP2];
  Nploc=3;
  Mploc=MT1; invMploc=invMT1;
endif
%
NdofP=Nploc*Nbtri;
%
NumPloc=linspace(1,Nploc*Nbtri,Nploc*Nbtri);
NumPloc=reshape(NumPloc,Nploc,Nbtri)';
%
Ku=sparse(Ndof,Ndof); Mu=sparse(Ndof,Ndof);
%
Mp=sparse(NdofP,NdofP); invMp=Mp;
%
Bx=sparse(NdofP,Ndof); By=Bx;
%
for t=1:Nbtri
   aire=Aires(t);
   AGLO=TriEdg(t,:);
   EdgNormT=EdgNorm(AGLO,:);
   Lga2T=Lga2(AGLO,:);
   for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
   end
   EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
   Lga2T(3)=EdgNormT(3,:)*EdgNormT(3,:)';
   [KNC,KCR,KFS]=KriNCTriangle(aire,EdgNormT,Lga2T,ordre);
   MuT=aire*MNC;
   %% Ku et Mu
   IGLO=NumUloc(t,:);
   Ku(IGLO,IGLO)+=KNC;
   Mu(IGLO,IGLO)+=MuT;
   % Matrice de couplage
   % -\int\div\vvec q ou \int\vvec\cdot\grad q
   [BNCx,BNCy,BCRx,BCRy,BFSx,BFSy]=BxyNCTriangle(aire,EdgNormT,ordre);
   PGLO=NumPloc(t,:);
   Bx(PGLO,IGLO)+=BNCx;
   By(PGLO,IGLO)+=BNCy;
   Mp(PGLO,PGLO)=aire*Mploc;
   invMp(PGLO,PGLO)=(1/aire)*invMploc;
end
