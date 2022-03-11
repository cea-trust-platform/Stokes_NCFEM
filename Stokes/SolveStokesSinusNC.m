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
% SolveStokes.m:
%
% Resout le probleme de Stokes "Sinus" 
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
% SYNOPSIS [Eu0,Eu1,Ep0,fig]=SolveStokesSinusNC
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - CoorNeu2(Nbpt+Nbedg,2) : coordonnees (x, y) des noeuds P2
%        - RefNeu(Nbpt,1) : reference des sommets
%        - CoorBary(Nbtri,3) :coordonnees (x, y) des barycentres des triangles
%        - CoorMil(Nbedg,2)   : Coordonnees des milieux d'aretes
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets) 
%        - NumTri2(4*Nbtri,3) : liste de triangles du maillage P2
%                   (3 numeros de sommets)
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - Lga2(Nbedg,1) : longueurs des aretes au carre
%        - EdgNorm(Nbedg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1) : aires des triangles
%        - fig : numero de figure si visualisation
%        - mi  : numéro du maillage
%        - ordre    : ordre d'approximation
% OUTPUT - Eu0 : Erreur L2 normalisee de la vitesse, calcul decompose
%        - Eu1 : Erreur H1 normaliseede la vitesse, calcul decompose
%        - Ep0 : Erreur L2 normalisee de la pression, calcul decompose
%        - fig = numero de la derniere figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Eu0,Eu1,Ep0,fig]=SolveStokesSinusNC(fig)
%
global Nbpt CoorNeu CoorNeu2 RefNeu
global Nbtri CoorBary Aires NumTri NumTri2 TriEdg
global Nbedg CoorMil RefEdg Lga2 EdgNorm EdgTri
global ordre mi RT algOct
%
[Ku,Mu,Mp,invMp,Bx,By] = MatStokesNC();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second membre  et solution exacte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (RT<=0)
  [RHS_Ux,RHS_Uy]=RHSStokesSinusNC();
else
  [RHS_Ux,RHS_Uy]=RHSStokesSinusRT();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Traitement des CL de Dirichlet
% Traitement des pressions constantes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DofU_A=find(RefEdg==0);
if (ordre==1)
  DofU=DofU_A;
  DofP=linspace(2,Nbtri,Nbtri-1)';
end
if (ordre==2)
  DofU_S=find(RefNeu==0);
  DofU_T=linspace(1,Nbtri,Nbtri)';
  NP2=Nbpt+Nbedg;
  DofU=[DofU_S;Nbpt+DofU_A;NP2+DofU_T];
  DofP=linspace(2,3*Nbtri,3*Nbtri-1)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolution du systeme lineaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndofP=size(Mp,1);
UnP=ones(ndofP,1);
Vol=UnP'*Mp*UnP;
RHS_p=zeros(ndofP,1);
[Uxex,Uyex,Pex,Uex2h,GUex2h,Pex2h]=SolutionStokesSinusNC(Mu,Ku,Mp);
if (algOct==1)
  [Uxh,Uyh,Ph]=DirectSolver(Ku,Bx,By,Mp,UnP,RHS_Ux,RHS_Uy,DofU,DofP,Vol);
else
  [Uxh,Uyh,Ph]=UzawaGCPMoy0(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,RHS_p,DofU,Vol);
end
%%%%%%%%%%%%%%%%%%%%
% Calcul des erreurs
%%%%%%%%%%%%%%%%%%%%
dUx=Uxh-Uxex; dUy=Uyh-Uyex; dP=Ph-Pex;
Ehu0=sqrt((dUx'*Mu*dUx+dUy'*Mu*dUy)/Uex2h);
Ehu1=sqrt((dUx'*Ku*dUx+dUy'*Ku*dUy)/GUex2h);
Ehp0=sqrt((dP'*Mp*dP)/Pex2h);
fprintf('P%i-NC mesh_%i, ||Ph(Uex)-U_h||_0/||Ph(Uex)||_0 = %7.2e\n',ordre,mi,Ehu0); 
fprintf('P%i-NC mesh_%i, ||Ph(Uex)-U_h||_h/||Ph(Uex)||_1 = %7.2e\n',ordre,mi,Ehu1);
fprintf('P%i-NC mesh_%i, ||Ph(Pex)-p_h||_0/||Ph(Pex)||_0 = %7.2e\n',ordre,mi,Ehp0); 
fprintf('------------------------------------------------\n');  
%%%%%%%%%%%%%%%%%%%%
% Calcul des erreurs
%%%%%%%%%%%%%%%%%%%%
[Eu0,Eu1,Ep0] = ErreursStokesSinusNC(Ku,Mu,Mp,Uxh,Uyh,Ph);
fprintf('P%i-NC mesh_%i, ||Uex-U_h||_0/||Uex||_0 = %7.2e\n',ordre,mi,Eu0); 
fprintf('P%i-NC mesh_%i, ||Uex-U_h||_h/||Uex||_1 = %7.2e\n',ordre,mi,Eu1);
fprintf('P%i-NC mesh_%i, ||Pex-P_h||_0/||pex||_0 = %7.2e\n',ordre,mi,Ep0);  
fprintf('------------------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fig>0)
  %% Matrice de masse LG
  [Uxh_LG,Uyh_LG,Ph_LG] = NCtoLG(Uxh,Uyh,Ph);
  %
  [Uxex_LG,Uyex_LG,Pex_LG] = SolutionStokesSinusLG(ordre);
  %
  if (ordre==1)
    NT=NumTri; CN=CoorNeu;
  end
  if (ordre==2)
    NT=NumTri2; CN=CoorNeu2; 
  end
  texX=sprintf('Solution exacte Ux, P%i, mesh%i', ordre,mi);
  tNCX=sprintf('Solution NC-FEM Ux, P%i, mesh%i', ordre,mi);
  %
  texY=sprintf('Solution exacte Uy, P%i, mesh%i', ordre,mi);
  tNCY=sprintf('Solution NC-FEM Uy, P%i, mesh%i', ordre,mi);
  %
  texP=sprintf('Solution exacte P, P%i, mesh%i', ordre,mi);
  tNCP=sprintf('Solution NC-FEM P, P%i, mesh%i', ordre,mi);
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Uxex_LG);
  view(2);
  shading interp
  title(texX)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),Uxh_LG);
  view(2);
  shading interp
  title(tNCX)
  colorbar;
  %
  fig=fig+1;
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Uyex_LG);
  view(2);
  shading interp
  title(texY)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),Uyh_LG);
  view(2);
  shading interp
  title(tNCY)
  colorbar;
  %
  fig=fig+1;
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Pex_LG);
  view(2);
  shading interp
  title(texP)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),Ph_LG);
  view(2);
  shading interp
  title(tNCP)
  colorbar;
end
