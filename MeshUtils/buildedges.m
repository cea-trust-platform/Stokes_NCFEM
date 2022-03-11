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
% buildedges.m:
% calcul des aretes et faces normales.
% 
% SYNOPSIS [CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,...
%          EdgTri,SomOpp,Lga2,EdgNorm,Aires]=buildedges(CoorNeu,RefNeu,NumTri,NumEdgB,RefEdgB)
%          
% INPUT  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - NumEdgB(NbEdgB,2) : Numero des 2 noeuds de chaque arete
%		     - RefEdgB(NbEdgB,1) : Reference de chaque arete
%
% OUTPUT - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des elements
%        - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete
%        - CoorMil(NbEdg,2)   : Coordonnees des milieux d'aretes
%		     - RefEdg(NbEdg,1) : Reference de chaque arete 
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
function [CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,...
          EdgTri,SomOpp,Lga2,EdgNorm,Aires]=buildedges(CoorNeu,RefNeu,NumTri,NumEdgB,RefEdgB)
% Lecture fichier maillage P1
% Algorithme pour trouver les aretes interieures
SomEdg = [1,2,3;2,3,1;3,1,2] ;
NumEdg= NumEdgB ;
NbEdgB=size(NumEdgB,1);
RefEdg= RefEdgB ;
NbEdg = NbEdgB  ; % nombre d'arete ;
Nbtri=size(NumTri,1);
% Pour chaque triangle, on donne le numero de l'arete opposee au sommet
TriEdg=zeros(Nbtri,3);
% Pour chaque arete, on donne le numero des triangles auquels elle appartient
EdgTri=zeros(NbEdg,2);
% Pour chaque arete, on donne le numero du sommet oppose a l'arete dans chaque triangle
%   SomOpp(a,2) = 0 si a est sur la frontiere 
SomOpp=zeros(NbEdg,2);
Aires=zeros(Nbtri,1);
CoorBary=zeros(Nbtri,2);
for t=1:Nbtri
  for i=1:3
    iglo=NumTri(t,i);
    CoorBary(t,:)=CoorBary(t,:)+CoorNeu(iglo,:);
  end
  CoorBary(t,:)=CoorBary(t,:)/3;
  % Vectors of edges 3 and 2
  S12=CoorNeu(NumTri(t,2),:)-CoorNeu(NumTri(t,1),:); % (X2-X1,Y2-Y1)
  S13=CoorNeu(NumTri(t,3),:)-CoorNeu(NumTri(t,1),:); % (X3-X1,Y3-Y1)
  % Jacobian of the transformation
  Jac=S12(1)*S13(2)-S12(2)*S13(1);
  % Volume of the element
  Aires(t)=0.5*abs(Jac);
  % Loop on the edges
  for aloc=1:3
    iloc = SomEdg(aloc,1);
    jloc = SomEdg(aloc,2);
    kloc = SomEdg(aloc,3);
    i=NumTri(t,iloc);   % un bout l'arete
    j=NumTri(t,jloc);   % autre bout de l'arete
    % numero des extremités par ordre croissant
    if (i>j)
      temp=i; i=j; j=temp;
    end
    k=NumTri(t,kloc);   % le numero de sommet oppose a l'arete
    xi=CoorNeu(i,1); yi=CoorNeu(i,2);
    xj=CoorNeu(j,1); yj=CoorNeu(j,2);
    exists=0;
    % Si et Sj sont des noeuds du bord
    if ( (RefNeu(i)~=0) && (RefNeu(j)~=0))
      % on cherche si SiSJ est une arete du bord
      for a=1:NbEdgB              
        if ((NumEdg(a,1)==i) && (NumEdg(a,2)==j))
          TriEdg(t,kloc)=a;
          EdgTri(a,1)=t;
          SomOpp(a,1)=k;
          exists=1;
          break;
        end                
      end
      if (exists==0)
        % Si l'arete SiSj n'existe pas comme arete du bord, 
        % SiSj devrait etre une arete interieure et
        % Sk devrait etre aussi un noeud du bord
        if (RefNeu(k)~=0)
          % on cherche si l'arete existe deja en tant qu'arete interieure
          if (NbEdg>NbEdgB)
            for a=NbEdgB+1:NbEdg
              if ((NumEdg(a,1)==i) && (NumEdg(a,2)==j))
                % l'arete existe deja
                TriEdg(t,kloc)=a;
                EdgTri(a,2)=t;
                SomOpp(a,2)=k;
                exists=1;
                break;
              end
            end
          end
        else
          printf('Que se passe-t-il ?\n');
        end
      end           
    else
      % SiSJ est une arete interieure
      for a=NbEdgB+1:NbEdg
        if ((NumEdg(a,1)==i) && (NumEdg(a,2)==j))
          TriEdg(t,kloc)=a;
          EdgTri(a,2)=t;
          SomOpp(a,2)=k;
          exists=1;
          break;
        end
      end
   end
   % on cree l'arete si elle n'existe pas
   if (exists==0)
     NbEdg=NbEdg+1;
     NumEdg= [NumEdg;i,j];
     TriEdg(t,kloc)=NbEdg;
     EdgTri=[EdgTri;t,0]; % EdgTri(a,1)=t;
     SomOpp=[SomOpp;k,0]; % SomOpp(a,1)=k;
     RefEdg(NbEdg)=0;
     exists=1;
   end %
  end % fin boucle sur les aretes   
end % fin boucle sur les triangles
NbEdg = size(NumEdg,1);
% Pour chaque arete, on donne les coordonnees du milieu de l'arete
CoorMil=zeros(NbEdg,2);
% vecteurs face-normale, orientes tri1->tri2
EdgNorm=zeros(NbEdg,2);
% longueurs des aretes au carre
Lga2=zeros(NbEdg,1);
% boucle sur les aretes
for a=1:NbEdg
  i=NumEdg(a,1); j=NumEdg(a,2); k=SomOpp(a,1);  
  Si=CoorNeu(i,:); Sj=CoorNeu(j,:);  Sk=CoorNeu(k,:);
  % (Si,Sj,Sk) dans le sens trigo ?
  SiSj=[Sj-Si,0]; SiSk=[Sk-Sj,0];
  z=cross(SiSj,SiSk);
  if (z(3)<0)
   NumEdg(a,1)=j;
   NumEdg(a,2)=i;
   i=NumEdg(a,1); j=NumEdg(a,2); 
   Si=CoorNeu(i,:); Sj=CoorNeu(j,:);
  end    
  EdgNormA=[Sj(2)-Si(2),Si(1)-Sj(1)];
  CoorMil(a,:)=0.5*(Si+Sj);
  SkMk=CoorMil(a,:)-CoorNeu(k,:);
  if (dot(EdgNormA,SkMk)<0)
    EdgNorm(a,:)=-EdgNormA;
  else
    EdgNorm(a,:)=EdgNormA;
  end
  Lga2(a)=dot(EdgNormA,EdgNormA);
end
for a=1:NbEdg
  EdgNormA=EdgNorm(a,:);
  Lga2(a,1)=dot(EdgNormA,EdgNormA);
end
% VERIFICATION
OK=0;
eps=1e-7;
if (OK)
for t=1:Nbtri
  X=zeros(3,3);
  for i=1:3
    X(i,1:2)=CoorNeu(NumTri(t,i),:);
  end
  % sens trigo ?
  %S1S2=[X(2,:)-X(1,:)];
  %S1S3=[X(3,:)-X(1,:)];
  %z=cross(S1S2,S1S3);
  %sens=1;
  %if (z(3)<0)
  %   printf('Triangle %i sens non trigo\n',t);
  %   sens=-1;
  %end
  %X(:,3)=ones(3,1);
  %Phi=-sens*2*Aires(t)*inv(X)'; 
  % GMSH : les sommets des triangles sont ranges dans le sens trigo
  X(:,3)=ones(3,1);   
  Phi=-2*Aires(t)*inv(X)';
  EdgNormT=TriEdg(t,:);
  for iloc=1:3
    ai = EdgNormT(iloc,:);
    if (EdgTri(ai,1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
    end
  end
  Diff=(EdgNormT-Phi(:,1:2));
  eDiff=sqrt(dot(Diff,Diff))/Lga(a);
  if (eDiff>eps)
    printf('Tri=%i,eDiff=%7.2e\n',t,eDiff);
  end
end
end
