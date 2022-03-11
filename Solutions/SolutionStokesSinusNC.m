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
% SolutionStokesSinusNC.m:
% Solution du probleme de Stokes "Sinus" 
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
% SYNOPSIS [Uxex,Uyx,Pex,Uex2h,GUex2h,Pex2h]=SolutionStokesSinusNC(Mu,Ku,Mp)
%
% GLOBAL - Nbtri            : nb de triangles
%        - CoorNeu(Nbpt,2)  : coordonnees (x, y) des sommets
%        - CoorMil(Nbedg,2) : coordonnees (x, y) des milieux d'aretes
%        - CoorBary(Nbtri,2): coordonnees (x, y) des barycentres des triangles
%        - NumTri(Nbtri,3)  : liste de triangles (3 numeros de sommets)
%        - Nbedg            : nb d'aretes
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
% OUTPUT - [Uxex,Uyex,Pex] : la solution exacte au points de discretisation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxex,Uyex,Pex,Uex2h,GUex2h,Pex2h]=SolutionStokesSinusNC(Mu,Ku,Mp)
global Nbpt CoorNeu RefNeu
global Nbtri CoorBary NumTri Aires
global Nbedg CoorMil TriEdg RefEdg
global ordre npi
% ordre1
if (ordre==1)
  Ndof_u=Nbedg;
  Ndof_p=Nbtri;
end
% ordre2
if (ordre==2)
  NP2=Nbpt+Nbedg;
  Ndof_u=NP2+Nbtri;
  Ndof_p=3*Nbtri;
end
Uxex=zeros(Ndof_u,1); Uyex=Uxex;
Pex=zeros(Ndof_p,1);
%
if (ordre==1)
  XY=npi*CoorMil;
  Uxex=(1-cos(XY(:,1))).*sin(XY(:,2));
  Uyex=(cos(XY(:,2))-1).*sin(XY(:,1));
  Pex =sin(npi*CoorBary(:,1)).*sin(npi*CoorBary(:,2));
end
%
if (ordre==2)
  PexLG=sin(npi*CoorNeu(:,1)).*sin(npi*CoorNeu(:,2));
  Pex=zeros(3*Nbtri,1);
  pglo=1;
  % Attention, c'est pas vraiment du P2
  ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
  MassUx=zeros(Ndof_u,1); MassUy=MassUx;
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    AGLO=Nbpt+TriEdg(t,:);
    TGLO=NP2+t;
    CoorNeuT=CoorNeu(IGLO,:);
    aire=Aires(t);
    pglof=pglo+2;
    Pex(pglo:pglof,1)=PexLG(IGLO,1);
    pglo=pglof+1;
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    npix=npi*xyp(:,1); npiy=npi*xyp(:,2);
    cosx=cos(npix); cosy=cos(npiy);
    sinx=sin(npix); siny=sin(npiy);
    awp=aire*wp';
    Uxloc=awp.*(1-cosx).*siny; 
    Uyloc=awp.*(cosy-1).*sinx;
    lambda2=lambda.*lambda;
    PhiS=2*lambda2-lambda;    
    PhiM(:,ILOC)=4*lambda(:,JLOC).*lambda(:,KLOC);
    PhiT=2-3*sum(lambda2,2);
    UGLO=[IGLO,AGLO,TGLO];
    PhiNC=[PhiS,PhiM,PhiT];
    MassUx(UGLO)+=PhiNC'*Uxloc;
    MassUy(UGLO)+=PhiNC'*Uyloc;
  end
  Uxex=Mu\MassUx; Uyex=Mu\MassUy;
end
% Traitement des CL de Dirichlet
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
GUex2h=Uxex'*Ku*Uxex+Uyex'*Ku*Uyex;
Uex2h=Uxex'*Mu*Uxex+Uyex'*Mu*Uyex; 
Pex2h=Pex'*Mp*Pex;
end
