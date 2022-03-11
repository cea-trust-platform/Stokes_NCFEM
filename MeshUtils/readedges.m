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
% readedges.m:
% routine de lecture des fichiers de maillages triangulaires 2D au format .edg 
% 
% SYNOPSIS [CoorBary,NumEdg,CoorMil,RefEdg,...
%             TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires]=readedges(filename)
%          
% INPUT  - filename : le nom d'un fichier au format edg sans son suffixe .edg
%
% OUTPUT - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des elements
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque Edg
%        - CoorMil(NbEdg,2)   : Coordonnees des milieux d'Edg
%		     - RefEdg(NbEdg,1) : Reference de chaque Edg 
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'Edg opposee au sommet Numtri(l,i)
%                  (3 numeros des Edg - matrice entiere Nbtri x 3)
%		     - EdgTri(NbEdg,2) : Pour chaque Edg, EdgTri(a,:) donne les numeros des 2 triangles de chaque Edg 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere 
%		     - SomOpp(NbEdg,2) : Numero du sommet oppose a l'Edg dans chaque triangle
%                                  SomOpp(a,2) = 0 si a est sur la frontiere 
%        - Lga2(NbEdg,1) : longueurs des Edg au carre
%        - EdgNorm(NbEdg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1) : aires des triangles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires]=readedges(filename)

edgesfile=strcat(filename,".edg");
fid2=fopen(edgesfile,'r');
if fid2 <=0,
  msg=['Le fichier de maillage : ' edgesfile ' n''a pas ete trouve'];
  error(msg);
end
while ~strcmp(fgetl(fid2),'$NbAretes'), end
NbEdg = str2num(fgetl(fid2));
% Pour chaque Edg, on donne les coordonnees du milieu de l'Edg
CoorMil=zeros(NbEdg,2);
% vecteurs face-normale, orientes tri1->tri2
EdgNorm=zeros(NbEdg,2);
% longueurs des Edg au carre
Lga2=zeros(NbEdg,1);
% boucle sur les Edg
for a=1:NbEdg
  tmp= str2num(fgetl(fid2));
  CoorMil(a,:) = tmp(2:3);
  EdgNorm(a,:) = tmp(4:5);
  Lga2(a)=tmp(6);
end
Lga=sqrt(Lga2);
while ~strcmp(fgetl(fid2),'$NumAretes'), end
% Pour chaque Edg, on donne les numeros de ses extremites
NumEdg=zeros(NbEdg,2);
% Pour chaque Edg, on indique si elle est au bord ou pas
RefEdg=zeros(NbEdg,1);
% Pour chaque Edg, on donne le numero des triangles auquels elle appartient
EdgTri=zeros(NbEdg,2);
% Pour chaque Edg, on donne le numero du sommet oppose a l'Edg dans chaque triangle
%   SomOpp(a,2) = 0 si a est sur la frontiere 
SomOpp=zeros(NbEdg,2);
% boucle sur les Edg
for a=1:NbEdg
  tmp= str2num(fgetl(fid2));
  RefEdg(a,1)=tmp(2);
  NumEdg(a,:)=tmp(3:4);
  EdgTri(a,:)=tmp(5:6);
  SomOpp(a,:)=tmp(7:8);
end
while ~strcmp(fgetl(fid2),'$Nbtri'), end
fprintf(fid2,'$Nbtri\n');
Nbtri = str2num(fgetl(fid2));
% Pour chaque triangle, on donne le numero de l'Edg opposee au sommet
TriEdg=zeros(Nbtri,3);
CoorBary=zeros(Nbtri,2);
Aires=zeros(Nbtri,1);
for t=1:Nbtri
  tmp= str2num(fgetl(fid2));
  TriEdg(t,:)=tmp(2:4);
  CoorBary(t,:)=tmp(5:6);
  Aires(t,1)=tmp(7);
end
fclose(fid2);
