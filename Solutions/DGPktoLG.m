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
% DGPktoLG.m:
%
% Projection DG base Pk vers LG P1 ou P2
%
% SYNOPSIS Vh_LG=DGPktoLG(Mu_LG,Mp,Vh,ordreT)
%          
% INPUT  - MV_LG : matrice de masse
%        - Vh(Ndof,1) : vecteur a  projeter
%        - ordreT(Nbtri,1) : ordre par triangle
%
% OUTPUT - Vh_LG : vecteur projete
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vh_LG=DGPktoLG(Mk_LG,Vh,k_DG,k_LG)
%
global Nbtri NumTri Nbpt Nbedg TriEdg Aires
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle;
nlocLG=3*k_LG;
% ORDRE 1
if (k_LG==1)
  MT=MT1; Ndof=Nbpt; NumTriM=NumTri;
end
% ORDRE 2
if (k_LG==2)
  MT=MT2; Ndof=Nbedg+Nbpt; NumTriM=[NumTri,Nbpt+TriEdg];
end
%
PVh=zeros(Ndof,1);
ideb=0;
for t=1:Nbtri
  kT=k_DG(t);
  nlocDG=3*kT;
  aire=Aires(t);
  if (nlocDG==0)
    nlocDG=1;
  end
  ifin=ideb+nlocDG;
  VhT=Vh(ideb+1:ifin,1);
  ideb=ifin;
  if (kT==k_LG)  
     MT_loc=aire*MT;
  end
  %
  if ((kT==1)&&(k_LG==2))
     MT_loc=aire*MT12';
  end
  %
  if ((kT==0)&&(k_LG==1))
     MT_loc=aire*MT01';
  end
  %
  if ((kT==0)&&(k_LG==2))
     MT_loc=aire*MT02';
  end
  %
  PVhT=MT_loc*VhT;
  %
  JGLO=NumTriM(t,:);
  PVh(JGLO)=PVh(JGLO)+PVhT;
end
Vh_LG=Mk_LG\PVh;
