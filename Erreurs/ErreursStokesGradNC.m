%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% ErreursStokesGradNC.m:
% Solution du probleme de Stokes  "Grad" 
%   avec les elements finis non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= 3[x^2,y^2]
%    div U = 0;
%
% U(x,y)=0; 
% p(x,y)=x^3+y^3-1/2;
%
% SYNOPSIS [Eu0,Eu1,Ep0] = ErreursStokesGradNC(Ku,Mu,Mp,Uxh,Uyh,Ph);
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3) : numeros des sommets des triangles
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2      
% INPUT  - Ku(Nu,Nu) matrice de raideur de la vitesse, composantes x ou y
%        - Mp : matrice de masse pression
%        - Uxh  : vitesse approchee, composante x
%        - Uyh  : vitesse approchee, composante x
%        - Ph : pression approchee
%
% OUTPUT - Eu0 : Erreur L2 vitesse normalisee, calcul decompose
%        - Eu1 : Erreur H1 vitesse normalisee, calcul decompose
%        - Ep0 : Erreur L2 pression normalisee, calcul decompose
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eu0,Eu1,Ep0] = ErreursStokesGradNC(Ku,Mu,Mp,Uxh,Uyh,Ph)
% 
global Nbpt CoorNeu 
global Nbtri Aires  NumTri 
global ordre
% Normes L2 de Uexact et grad(Uexact)
Pex2=9/56;
%
Tri=linspace(1,Nbtri,Nbtri); Tri=Tri';
%
if (ordre==1)
  NumPloc=Tri;
endif
%
if (ordre==2)
  NumPloc=linspace(1,3*Nbtri,3*Nbtri);
  NumPloc=reshape(NumPloc,3,Nbtri)';
endif
%
% Erreur ||Pex-Ph||_0
Ep0=0;
PhPex=0;
PhiM=zeros(7,3);
for t=1:Nbtri
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  %
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  %
  Pex=awp.*(xyp(:,1).^3+xyp(:,2).^3-0.5);
  %
  % CALCULS DE (Uh,Uex)_{0,T} (Ph,pex)_{0,T} (\grad Uh,\grad Uex)_{0,T}
  if (ordre==1)
    PhPex+=Ph(t,1)*sum(Pex);
  end
  if (ordre==2)
    PGLO=NumPloc(t,:); PT=Ph(PGLO,1);
    PhPex+=Pex'*(lambda*PT);
  end
end
% 
Eu0=sqrt(Uxh'*Mu*Uxh+Uyh'*Mu*Uyh);
Eu1=sqrt(Uxh'*Ku*Uxh+Uyh'*Ku*Uyh);
Ep0=sqrt((Ph'*Mp*Ph-2*PhPex+Pex2)/Pex2);
%
