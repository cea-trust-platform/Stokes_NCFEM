%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% readP2mesh.m:
% routine de lecture des fichiers de maillages triangulaires 2D au format .P2 
% 
% SYNOPSIS [CoorNeu2,RefNeu2,NumTri2,RefTri2]=readP2mesh(filename)
%          
% INPUT  - filename : le nom d'un fichier au format P2 sans son suffixe .P2
%
% OUTPUT - CoorNeu2(Nbpt2,2) : coordonnees (x, y) des sommets P2
%        - RefNeu(Nbpt2,1) : reference des sommets
%        - NumTri(Nbtri2,3)  : liste de triangles P1 (3 numeros de sommets)
%        - NumTri2(Nbtri2,3) : liste de triangles P2 (3 numeros de sommets)
%        - RefTri2(Nbtri2,1) : Reference de chaque triangle maillage P2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoorNeu2,RefNeu2,NumTri2,RefTri2]=readP2mesh(filename)
%%%%%%%%%%%%%%%%%%
% P2
%%%%%%%%%%%%%%%%%%
P2file=strcat(filename,".P2");
fid3=fopen(P2file,'r');
if fid3 <=0,
  msg=['Le fichier de maillage : ' P2file ' n''a pas ete trouve'];
  error(msg);
end
while ~strcmp(fgetl(fid3),'$Nodes'), end
Nbpt2 = str2num(fgetl(fid3));
for i=1:Nbpt2
  tmp= str2num(fgetl(fid3));
  CoorNeu2(i,:)=tmp(2:3);
  RefNeu2(i,1)=tmp(4);
end
while ~strcmp(fgetl(fid3),'$Elements'), end
Nbtri2 = str2num(fgetl(fid3));
for t=1:Nbtri2
  tmp= str2num(fgetl(fid3));
  RefTri2(t,1)=tmp(2);
  NumTri2(t,:)=tmp(3:5);    
end
fclose(fid3);
