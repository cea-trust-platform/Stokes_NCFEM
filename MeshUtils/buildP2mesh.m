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
% buildP2mesh.m:
% Raffinement d'un maillage P1 avec 4 sous-triangles par triangles en utilisant les milieux d'aretes
% Utile pour la visu du P2 ou pour raffiner ulterieurement
%
% SYNOPSIS [CoorNeu2,RefNeu2,NumTri2,RefTri2]=buildP2mesh(CoorNeu,RefNeu,NumTri,RefTri,NumEdg,CoorMil,RefEdg,TriEdg
%          
% INPUT  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles du maillage P1
%        - RefTri(Nbtri,1) : reference des triangles du maillage P1
%        - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete du maillage P1
%        - CoorMil(NbEdg,2)   : Coordonnees des milieux d'aretes du maillage P1
%		     - RefEdg(NbEdg,1) : Reference de chaque arete du maillage P1
%		     - TriEdg(Nbtri,3) : Pour chaque triangle P1, Triarete(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%
% OUTPUT - CoorNeu2(Nbpt2,2) : coordonnees (x, y) des sommets suivis des milieux d'aretes
%        - RefNeu2(Nbpt2,1) : reference des sommets suivis des milieux d'aretes
%        - NumTri2(Nbtri2,3) : liste de triangles du maillage raffine (3 numeros de noeuds)
%        - RefTri2(Nbtri2,1) : reference des triangles du maillage raffine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoorNeu2,RefNeu2,NumTri2,RefTri2]=buildP2mesh(CoorNeu,RefNeu,NumTri,RefTri,NumEdg,CoorMil,RefEdg,TriEdg)
% 
Nbtri=size(NumTri,1);
Nbtri2=4*Nbtri;
Nbpt=size(CoorNeu,1);
NbEdg=size(NumEdg,1);
Nbpt2=Nbpt+NbEdg;

NumTri2=[];
RefTri2=[];

CoorNeu2=zeros(Nbpt2,2);
CoorNeu2(1:Nbpt,:)=CoorNeu;

RefNeu2=zeros(Nbpt2,1);
RefNeu2(1:Nbpt)=RefNeu;

SomEdg=[3,2;1,3;2,1];

triloc=zeros(1,3);
for t=1:Nbtri
  for aloc=1:3
    a=TriEdg(t,aloc);
    CoorNeu2(a+Nbpt,:)=CoorMil(a,:);
    RefNeu2(a+Nbpt)=0;
    % arete du bord ?
    if (RefEdg(a)~=0)       
       RefNeu2(a+Nbpt)=1;
    end
    triloc(aloc)=a+Nbpt;
  end
  NumTri2=[NumTri2;triloc];
  RefTri2=[RefTri2;RefTri(t)];
  for iloc=1:3
    jloc=SomEdg(iloc,1);
    kloc=SomEdg(iloc,2);
    s1=NumTri(t,iloc);
    m3=TriEdg(t,jloc)+Nbpt;
    m2=TriEdg(t,kloc)+Nbpt;
    triloc=[s1,m3,m2];
    NumTri2=[NumTri2;triloc];
    RefTri2=[RefTri2;RefTri(t)];
  end
end
NbtriP2=size(NumTri2,1);
if (Nbtri2~=NbtriP2)
   printf('Pb nombre de triangles, Nbtri2=%i, NbtriP2= %i\n',Nbtri2,NbtriP2);
end
% VERIFICATIONS
for t=1:Nbtri2
  X=zeros(3,3);
  for iloc=1:3
    i=NumTri2(t,iloc);
    X(iloc,1:2)=CoorNeu2(i,:);
  end
  % sens trigo ?
  S1S2=[X(2,:)-X(1,:)];
  S1S3=[X(3,:)-X(1,:)];
  z=cross(S1S2,S1S3);
  sens=1;
  if (z(3)<0)
   printf('Triangle P2 %i sens non trigo\n',t);
   sens=-1;
   % on echange les deux derniers elements
   temp=NumTri2(t,2);
   NumTri(t,2)=NumTri(t,3);
   NumTri(t,3)=temp;
  end
end
