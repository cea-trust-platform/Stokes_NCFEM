%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
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