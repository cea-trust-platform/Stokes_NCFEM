%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% SolveStokes.m:
%
% Resout le probleme de Stokes "Grad" 
%   avec les elements finis non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= 3[x^2,y^2]
%    div U = 0;
%
% U(x,y)=0
% p(x,y)=x^3+y^3-1/2;
%
% SYNOPSIS [Eu0,Eu1,Ep0,fig]=SolveStokesGradNC
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

function [Eu0,Eu1,Ep0,fig]=SolveStokesGradNC(fig)
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
if (RT==0)
  [RHS_Ux,RHS_Uy]=RHSStokesGradNC();
else
  [RHS_Ux,RHS_Uy]=RHSStokesGradRT();
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
[Pex,Pex2h]=SolutionStokesGradNC(Mp);
if (algOct==1)
  [Uxh,Uyh,Ph]=DirectSolver(Ku,Bx,By,Mp,UnP,RHS_Ux,RHS_Uy,DofU,DofP,Vol);
else
  [Uxh,Uyh,Ph]=UzawaGCPMoy0(Ku,Bx,By,Mp,invMp,UnP,RHS_Ux,RHS_Uy,RHS_p,DofU,Vol);
end
%%%%%%%%%%%%%%%%%%%%
% Calcul des erreurs
%%%%%%%%%%%%%%%%%%%%
[Eu0,Eu1,Ep0] = ErreursStokesGradNC(Ku,Mu,Mp,Uxh,Uyh,Ph);
dP=Ph-Pex; Ehp0=sqrt((dP'*Mp*dP)/Pex2h);
fprintf('P%i-NC mesh_%i, ||U_h||_0 = %7.2e\n',ordre,mi,Eu0); 
fprintf('P%i-NC mesh_%i, ||U_h||_h = %7.2e\n',ordre,mi,Eu1);
fprintf('P%i-NC mesh_%i, ||Pex-P_h||_0/||Pex||_0 = %7.2e\n',ordre,mi,Ep0);
fprintf('P%i-NC mesh_%i, ||Ph(Pex)-p_h||_0/||Ph(Pex)||_0 = %7.2e\n',ordre,mi,Ehp0);   
fprintf('------------------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fig>0)
  %% Matrice de masse LG
  [Uxh_LG,Uyh_LG,Ph_LG] = NCtoLG(Uxh,Uyh,Ph);
  %
  if (ordre==1)
    NT=NumTri; CN=CoorNeu;
  end
  if (ordre==2)
    NT=NumTri2; CN=CoorNeu2; 
  end
  %
  Pex_LG =CN(:,1).^3+CN(:,2).^3-0.5;
  %
  tNCX=sprintf('Solution NC-FEM Ux, P%i, mesh%i', ordre,mi);
  tNCY=sprintf('Solution NC-FEM Uy, P%i, mesh%i', ordre,mi);
  %
  texP=sprintf('Solution exacte P, P%i, mesh%i', ordre,mi);
  tNCP=sprintf('Solution NC-FEM P, P%i, mesh%i', ordre,mi);
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Uxh_LG);
  view(2);
  shading interp
  title(tNCX)
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