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
% MatMassNCtoLG.m:
% Matrice de projection NC vers LG P1 ou P2
%
% SYNOPSIS MatP=MatMassNCtoLG(Nbpt,NbEdg,NumTri,TriEdg,Aires,ordre)
%          
% INPUT  - Nbpt : nombre de sommets
%        - NbEdg : nombre d'aretes
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'Edg opposee au sommet NumTri(l,i)
%                  (3 numeros des Edg - matrice entiere Nbtri x 3)
%		     - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2 
%
% OUTPUT - MatP : Matrice de projection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatP=MatMassNCtoLG(Nbpt,NbEdg,NumTri,TriEdg,Aires,ordre)

Nbtri=size(NumTri,1);
%
if (ordre==1)
  MatP=sparse(Nbpt,NbEdg);
  for t=1:Nbtri
    aire6=Aires(t)/6;
    for iloc=1:3
      iglo=NumTri(t,iloc);
      for aloc=1:3
        if (aloc~=iloc)
           aglo=TriEdg(t,aloc);
           MatP(iglo,aglo)+=aire6;
        endif       
      endfor
    endfor
  end
endif

if (ordre==2)
  Ndof=Nbpt+NbEdg;
  MatP=[Mu_LG,sparse(Ndof,Nbtri)];
  TGLO=Ndof;
  for t=1:Nbtri
    aire5=Aires(t)/5;
    aire30=Aires(t)/30;
    TGLO=TGLO+1;
    IGLO=NumTri(t,:);
    MatP(IGLO,TGLO)=aire30*ones(3,1);
    AGLO=TriEdg(t,:)+Nbpt;
    MatP(AGLO,TGLO)=aire5*ones(3,1);
  end
endif
