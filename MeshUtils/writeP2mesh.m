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
