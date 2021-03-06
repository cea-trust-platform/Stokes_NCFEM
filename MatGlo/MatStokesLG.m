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
% MatStokesLG.m:
% Assemblages des matrices de masse, de raideur et de couplage 
%   EF de Taylor-Hood ou de Scott-Vogelius
%
% SYNOPSIS [Ku,Mu,Mp,invMp,Bx,By] = MatStokesLG(CoorNeu,NumTri,TriEdg,EdgTri,Lga2,EdgNorm,Aires,type,ordre,nuT)
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
%        - type : 'SV' pour Scott-Vogelius (P2-P1disc) ou 'TH' pour Taylor-Hood (P2-P1)
%        - nuT : coefficient de raideur
%
% OUTPUT - Ku : matrice de raideur vitesse
%        - Mu : matrice de masse vitesse
%        - Mp : matrice de masse pression
%        - invMp : matrice inverse de la matrice de masse pression
%        - Bx : matrice de couplage vitesse-pression
%        - By : matrice de couplage vitesse-pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu,Mp,invMp,Bx,By] = MatStokesLG(CoorNeu,NumTri,TriEdg,EdgTri,Lga2,EdgNorm,Aires,type,nuT)
%
ordre=2;
%
[Ku,Mu]=MatPoissonLG(CoorNeu,NumTri,TriEdg,EdgTri,Lga2,EdgNorm,Aires,ordre,nuT);
%
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle;
%
Nbpt =size(CoorNeu,1);
Nbtri=size(NumTri,1);
NbEdg=size(EdgTri,1);
NumUloc=[NumTri,Nbpt+TriEdg];
Ndof=Nbpt+NbEdg;
%-------------%
if (type=='TH')
  NumPloc=NumTri;
  NdofP=Nbpt;invMp=sparse(NdofP,1);
end
%-------------%
if (type=='SV')
  %
  NumPloc=linspace(1,3*Nbtri,3*Nbtri);
  NumPloc=reshape(NumPloc,3,Nbtri)'; 
  NdofP=3*Nbtri;  invMp=sparse(NdofP,NdofP);
end
Mp=sparse(NdofP,NdofP);
%
Bx=sparse(NdofP,Ndof); By=Bx;
%
for t=1:Nbtri
   aire=Aires(t);
   AGLO=TriEdg(t,:);
   EdgNormT=EdgNorm(AGLO,:);
   for iloc=1:2
       if (EdgTri(AGLO(iloc),1)~=t)
          EdgNormT(iloc,:)=-EdgNormT(iloc,:);
       end
   end
   EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
   % Matrice de couplage
   [B1x,B1y,B2x,B2y]=BxyLGTriangle(aire,EdgNormT,ordre);
   PGLO=NumPloc(t,:);
   IGLO=NumUloc(t,:);
   Bx(PGLO,IGLO)+=B2x;
   By(PGLO,IGLO)+=B2y;
   Mp(PGLO,PGLO)+=aire*MT1;
   if (type=='TH')
      invMp(PGLO)=invMp(PGLO)+aire;
   end
   %
   if (type=='SV')
     invMp(PGLO,PGLO)=(1/aire)*invMT1;
   end
end
%
if (type=='TH')
  % calcul de invMP=matrice de masse P1 reduite
  invMp=3*spdiags(1./invMp,0,Nbpt,Nbpt);
end
