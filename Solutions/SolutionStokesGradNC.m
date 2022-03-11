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
% SolutionStokesGradNC.m:
% Solution du probleme de Stokes  "Grad" 
%   avec les elements finis non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= 3[x^2,y^2]
%    div U = 0;
%
% U(x,y)=0
% p(x,y)=x^3+y^3-1/2;
%
% SYNOPSIS [Pex,Pex2h]=SolutionStokesGradNC(Mp)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - CoorMil(Nbedg,2) : coordonnees (x, y) des milieux d'aretes
%        - CoorBary(Nbtri,2): coordonnees (x, y) des barycentres des triangles
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Nbedg            : nb d'aretes
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
% OUTPUT - [Pex] : la solution exacte au points de discretisation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pex,Pex2h]=SolutionStokesGradNC(Mp)
global Nbpt CoorNeu
global Nbtri NumTri
global Nbedg CoorBary 
global ordre
% ordre1
if (ordre==1)
  Ndof_u=Nbedg;
  Ndof_p=Nbtri;
end
% ordre2
if (ordre==2)
  Ndof_u=Nbpt+Nbedg+Nbtri;
  Ndof_p=3*Nbtri;
end
Uxex=zeros(Ndof_u,1); Uyex=Uxex;
Uex2h=0; GUex2h=0;
Pex=zeros(Ndof_p,1);
%
if (ordre==1)
  Pex =CoorBary(:,1).^3+CoorBary(:,2).^3-0.5;
end
%
if (ordre==2)
  PexLG=CoorNeu(:,1).^3+CoorNeu(:,2).^3-0.5;
  Pex=zeros(3*Nbtri,1);
  pglo=1;
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    CoorNeuT=CoorNeu(IGLO,:);
    pglof=pglo+2;
    Pex(pglo:pglof,1)=PexLG(IGLO,1);
    pglo=pglof+1;
  end
end
Pex2h=Pex'*Mp*Pex;
end
