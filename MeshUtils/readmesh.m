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
% readmesh.m:
% routine de lecture de fichiers de maillages triangulaires 2D au format .msh
% et calcul des aretes et faces normales.
% 
% SYNOPSIS [CoorNeu,RefNeu,NumTri,RefTri,NumEdgB,RefEdgB]=readmesh(filename)
%          
% INPUT  - filename : le nom d'un fichier de maillage au format msh
%                   SANS SON SUFFIXE .msh
%
% OUTPUT - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - RefTri(Nbtri,1) : Reference de chaque triangle
%        - NumEdgB(NbEdgB,2) : Numero des 2 noeuds de chaque arete
%		     - RefEdgB(NbEdgB,1) : Reference de chaque arete 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoorNeu,RefNeu,NumTri,RefTri,NumEdgB,RefEdgB]=readmesh(filename)
% MAILLAGE INITIAL
meshfile=strcat(filename,'.msh');
fid=fopen(meshfile,'r');
if fid <=0,
  msg=['Le fichier de maillage : ' meshfile ' n''a pas ete trouve'];
  error(msg);
end
while ~strcmp(fgetl(fid),'$Nodes'), end
Nbpt = str2num(fgetl(fid));
CoorNeu = zeros(Nbpt,2);
RefNeu = zeros(Nbpt,1);
NbEdgB = 0;
% Numero des 2 noeuds de chaque arete
NumEdgB = [];
RefEdgB = [];
RefNeuBis = zeros(Nbpt,1);
for i=1:Nbpt
  tmp= str2num(fgetl(fid));
  CoorNeu(i,:) = tmp(2:3);
end
while ~strcmp(fgetl(fid),'$Elements'), end
Nbtri = str2num(fgetl(fid));
tmp= str2num(fgetl(fid)); 
test = tmp(2);
% Aretes et noeuds du bord
NbAretesBord=0;
while test==1
  RefNeuBis(tmp(6:7)')=tmp(2);
  i=tmp(6); j=tmp(7);
  if (i>j)
   temp=i; i=j; j=temp;
  end
  NumEdgB= [NumEdgB;i,j];
  RefEdgB= [RefEdgB;tmp(5)];    
  tmp= str2num(fgetl(fid));
  test = tmp(2);
  Nbtri = Nbtri-1;
  NbEdgB=NbEdgB+1;
end
RefNeu(find(RefNeu==0))=RefNeuBis(find(RefNeu==0));
% Triangles
NumTri = zeros(Nbtri,3);
RefTri = zeros(Nbtri,1);
for i=1:Nbtri
  NumTri(i,:) = tmp(end-2:end);
  RefTri(i)=tmp(4);
  % Milieu homogene
  % triangle du bord  
  tmp= str2num(fgetl(fid));
end
fclose(fid);
