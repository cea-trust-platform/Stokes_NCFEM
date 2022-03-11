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
% writeedges.m:
% routine d'ecriture de fichiers de maillages triangulaires 2D au format .edg 
% 
% SYNOPSIS writeedges(CoorMil,EdgNorm,Lga2,RefEdg,NumEdg,EdgTri,SomOpp,TriEdg,CoorBary,Aires,filename)
%          
% INPUT  - CoorBary(Nbtri,2) : coordonnees (x, y) des barycentres des elements
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
%        - filename : le nom d'un fichier de maillage SANS SUFFIXE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeedges(CoorBary,NumEdg,CoorMil,RefEdg,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires,filename)
%
% Enregistrement des donnees des Edg.
edgesfile=strcat(filename,".edg");
fid2=fopen(edgesfile,'w');
%
NbEdg=size(EdgNorm,1);
Nbtri=size(TriEdg,1);
%
fprintf(fid2,'$NbAretes\n');
fprintf(fid2,'%i\n',NbEdg);
% OK pour maillage d'un carre [0,1]x[0,1].
for a=1:NbEdg
    fprintf(fid2,'%i %2.18e %2.18e %2.18e %2.18e %2.18e\n',a,CoorMil(a,:),EdgNorm(a,:),Lga2(a,1));
end
fprintf(fid2,'$EndNbAretes\n');
fprintf(fid2,'$NumAretes\n');
for a=1:NbEdg
    fprintf(fid2,'%i %i %i %i %i %i %i %i\n',a,RefEdg(a,1),NumEdg(a,:),EdgTri(a,:),SomOpp(a,:));
end
fprintf(fid2,'$EndNumAretes\n');
fprintf(fid2,'$Nbtri\n');
fprintf(fid2,'%i\n',Nbtri);
for t=1:Nbtri
    fprintf(fid2,'%i %i %i %i %2.18e %2.18e %2.18e\n',t,TriEdg(t,:),CoorBary(t,:),Aires(t,1));
end
fprintf(fid2,'$EndNbtri\n');
fclose(fid2);
