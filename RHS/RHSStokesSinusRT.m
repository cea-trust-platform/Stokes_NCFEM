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
% RHSStokesSinusRT.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "Sinus" 
%   avec les elements finis non conformes CR ou FS 

%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= npi*[ sin(npi*y)*{(b-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                                sin(npi*x)*{(b+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=sin(npi*x)*sin(npi*y);
%
%
% SYNOPSIS [RHS_Ux,RHS_Uy]=RHSStokesSinusRT(Mu)
%          
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles (3 numeros de sommets)
%        - Nbedg : nombre d'aretes
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - RefEdg(Nbedg,1) : Reference de chaque arete
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2
%        - nu : viscosite
%
% OUTPUT - RHS_Ux : second membre composante x
%        - RHS_Uy : second membre composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesSinusRT()
%
global Nbpt CoorNeu RefNeu
global Nbtri Aires NumTri
global Nbedg TriEdg RefEdg EdgTri EdgNorm
global nu ordre
global npi npi2
%
npi2nu=npi2*nu;
cx=npi-2*npi2nu;
cy=npi+2*npi2nu;
Ndof=Nbedg;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  Ndof=NP2+Nbtri;
end
RHS_Ux=zeros(Ndof,1); RHS_Uy=zeros(Ndof,1);
%
% EF de CR P1-NC (ordre=1) 
if (ordre==1)
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    CoorNeuT=CoorNeu(IGLO,:);
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    awp=0.5*wp';
    sinxy=sin(npi*xyp); cosxy=cos(npi*xyp);
    FxLOC = awp.*sinxy(:,2).*(cx*cosxy(:,1)+npi2nu); 
    FyLOC = awp.*sinxy(:,1).*(cy*cosxy(:,2)-npi2nu);
    % (d|T|)^{-1}\int_T(x-OSi)\cdot\Fvec
    RHS_tmp=(ones(3,1)*xyp(:,1)'-CoorNeuT(:,1)*ones(1,np))*FxLOC+...
            (ones(3,1)*xyp(:,2)'-CoorNeuT(:,2)*ones(1,np))*FyLOC;
    %
    % vecteurs face-normale
    AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
    for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
    RHS_Ux(AGLO)+=EdgNormT(:,1).*RHS_tmp;
    RHS_Uy(AGLO)+=EdgNormT(:,2).*RHS_tmp;
    %
  end
end
%
% EF de FS (ordre=2) 
if (ordre==2)
  ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
  phiM=zeros(7,3);
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    CoorNeuT =CoorNeu(IGLO,:);
    % vecteurs face-normale
    AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
    for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
    % volume of the element
    aire=Aires(t);
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    awp=aire*wp';
    sinxy=sin(npi*xyp); cosxy=cos(npi*xyp);
    FxLOC = awp.*sinxy(:,2).*(cx*cosxy(:,1)+npi2nu); 
    FyLOC = awp.*sinxy(:,1).*(cy*cosxy(:,2)-npi2nu);
    FLOC = [FxLOC,FyLOC];
    %
    [PiRTxx,PiRTxy,PiRTyx,PiRTyy]=PiRT1(aire,CoorNeuT,EdgNormT,lambda);
    %
    UGLO=[IGLO,AGLO+Nbpt,NP2+t];
    RHS_Ux(UGLO)+=(FxLOC'*PiRTxx+FyLOC'*PiRTxy)';
    RHS_Uy(UGLO)+=(FxLOC'*PiRTyx+FyLOC'*PiRTyy)';
    % 
  end
end
%
%% Traitement des CL de Dirichlet
Dof_DiriA=find(RefEdg==0);
Dof_BordA=find(RefEdg~=0);
if (ordre==1)
  Dof_Diri=Dof_DiriA;
  Dof_Bord=Dof_BordA;
end
if (ordre==2)
  Dof_DiriS=find(RefNeu==0);
  Dof_BordS=find(RefNeu~=0);
  Dof_DiriT=linspace(1,Nbtri,Nbtri)';
  NP2=Nbpt+Nbedg;
  Dof_Diri=[Dof_DiriS;Nbpt+Dof_DiriA;...
            NP2+Dof_DiriT];
  Dof_Bord=[Dof_BordS;Nbpt+Dof_BordA];
end
Nbord=size(Dof_Bord,1);
RHS_Ux(Dof_Bord,1)=zeros(Nbord,1);
RHS_Uy(Dof_Bord,1)=zeros(Nbord,1);
%
