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
% readmeshfiles.m:
% routine de lecture des fichiers de maillages triangulaires 2D au format .msh 
% 
% SYNOPSIS [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,CoorMil,...
%             RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires]=readmeshfiles(filename)
%          
% INPUT  - filename : le nom d'un fichier sans suffixe
%
% OUTPUT - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets P1
%        - CoorNeu2(Nbpt,2) : coordonnees (x, y) des sommets P2
%        - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des elements
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3)  : liste de triangles P1 (3 numeros de sommets)
%        - NumTri2(Nbtri,3) : liste de triangles P2 (3 numeros de sommets)
%        - RefTri(Nbtri,1)  : Reference de chaque triangle maillage P1
%        - RefTri2(Nbtri,1) : Reference de chaque triangle maillage P2
%        - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete
%        - NumEdgB(NbEdgB,2) : Numero des 2 noeuds de chaque arete du bord
%        - CoorMil(NbEdg,2)   : Coordonnees des milieux d'aretes
%		     - RefEdg(NbEdg,1) : Reference de chaque arete 
%		     - RefEdgB(NbEdgB,1) : Reference de chaque arete du bord
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(NbEdg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere 
%		     - SomOpp(NbEdg,2) : Numero du sommet oppose a l'arete dans chaque triangle
%                                  SomOpp(a,2) = 0 si a est sur la frontiere 
%        - Lga2(NbEdg,1) : longueurs des aretes au carre
%        - EdgNorm(NbEdg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1) : aires des triangles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,NumEdgB,CoorMil,...
           RefEdg,RefEdgB,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires]=readmeshfiles(filename)
%%%%%%%%%%%%%%%%%%
% MAILLAGE INITIAL
%%%%%%%%%%%%%%%%%%
[CoorNeu,RefNeu,NumTri,RefTri,NumEdgB,RefEdgB]=readmesh(filename);
%%%%%%%%%%%%%%%%%%
% ARETES
%%%%%%%%%%%%%%%%%%
[CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires]=readedges(filename);
%%%%%%%%%%%%%%%%%%
% P2
%%%%%%%%%%%%%%%%%%
[CoorNeu2,RefNeu2,NumTri2,RefTri2]=readP2mesh(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coherence vecteurs face-normale au bord 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NbEdg=size(RefEdg,1);
for f=1:NbEdg
  if (RefEdg(f)~=0)
    t=EdgTri(f,1);
    AGLO=TriEdg(t,:);
    EdgNormT=EdgNorm(AGLO,:);
    numloc=3;
    for floc=1:2
       if (AGLO(floc)==f)
         numloc=floc;
       end
       %
       if (EdgTri(AGLO(floc),1)~=t)
         EdgNormT(floc,:)=-EdgNormT(floc,:);
       end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);    
    EdgNorm(f,:)=EdgNormT(numloc,:);
  end
end
