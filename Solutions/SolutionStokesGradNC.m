%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
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