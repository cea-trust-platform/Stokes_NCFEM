%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% writeP2mesh.m:
% routine d'ecriture de fichiers de maillages triangulaires 2D au format .P2 
% 
% SYNOPSIS writeP2mesh(CoorNeu2,RefNeu2,NumTri2,RefTri2,filename)
%          
% INPUT  - CoorNeu2(Nbpt2,2) : coordonnees (x, y) des sommets suivis des milieux d'aretes
%        - RefNeu2(Nbpt2,1) : reference des sommets suivis des milieux d'aretes
%        - NumTri2(Nbtri2,3) : liste de triangles du maillage raffine (3 numeros de noeuds)
%        - RefTri2(Nbtri2,1) : reference des triangles du maillage raffine
%        - filename : le nom d'un fichier de maillage SANS SUFFIXE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeP2mesh(CoorNeu2,RefNeu2,NumTri2,RefTri2,filename)
%
Nbpt2=size(CoorNeu2,1);
Nbtri2=size(NumTri2,1);
%
P2file=strcat(filename,".P2");
fid=fopen(P2file,'w');
% Enregistrement du maillage P2
fprintf(fid,'$Nodes\n');
fprintf(fid,'%i\n',Nbpt2);
for i=1:Nbpt2
  fprintf(fid,'%i %2.18e %2.18e %i\n',i,CoorNeu2(i,:),RefNeu2(i,1));
end
fprintf(fid,'$EndNodes\n');
fprintf(fid,'$Elements\n');
fprintf(fid,'%i\n',Nbtri2);
for t=1:Nbtri2
  fprintf(fid,'%i %i %i %i %i\n',t,RefTri2(t,1),NumTri2(t,:));
end
fprintf(fid,'$EndElements\n');
fclose(fid);
end