%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% MatStokesNC.m:
%
% Assemblages des matrices de masse, de raideur et de couplage
%             EF de CR P1-NC (ordre=1) ou FS (ordre=2) 
%
% SYNOPSIS [Ku,Mu,Mp,invMp,Bx,By] = MatStokesNC
%          
% INPUT  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - Lga2(Nbedg,1)    : longueurs des aretes au carre
%        - EdgNorm(Nbedg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - Ku : matrice de raideur vitesse
%        - Mu : matrice de masse vitesse
%        - Mp : matrice de masse pression
%        - invMp : inverse de la matrice de masse pression
%        - Bx : matrice de couplage vitesse-pression
%        - By : matrice de couplage vitesse-pression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu,Mp,invMp,Bx,By] = MatStokesNC()
%
global Nbpt CoorNeu
global Nbtri Aires NumTri TriEdg
global Nbedg Lga2 EdgNorm EdgTri
global ordre
%
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle();
[MNC,MCR,MFS]=MassNCTriangle(ordre,MT2);
if (ordre==1)
  Ndof=Nbedg;
  NumUloc=TriEdg;
  Nploc=1;
  Mploc=1; invMploc=1;
endif
%
if (ordre==2)
  NP2=Nbpt+Nbedg;
  Ndof=NP2+Nbtri;
  Tri=linspace(1,Nbtri,Nbtri)';
  NumUloc=[NumTri,Nbpt+TriEdg,Tri+NP2];
  Nploc=3;
  Mploc=MT1; invMploc=invMT1;
endif
%
NdofP=Nploc*Nbtri;
%
NumPloc=linspace(1,Nploc*Nbtri,Nploc*Nbtri);
NumPloc=reshape(NumPloc,Nploc,Nbtri)';
%
Ku=sparse(Ndof,Ndof); Mu=sparse(Ndof,Ndof);
%
Mp=sparse(NdofP,NdofP); invMp=Mp;
%
Bx=sparse(NdofP,Ndof); By=Bx;
%
for t=1:Nbtri
   aire=Aires(t);
   AGLO=TriEdg(t,:);
   EdgNormT=EdgNorm(AGLO,:);
   Lga2T=Lga2(AGLO,:);
   for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
   end
   EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
   Lga2T(3)=EdgNormT(3,:)*EdgNormT(3,:)';
   [KNC,KCR,KFS]=KriNCTriangle(aire,EdgNormT,Lga2T,ordre);
   MuT=aire*MNC;
   %% Ku et Mu
   IGLO=NumUloc(t,:);
   Ku(IGLO,IGLO)+=KNC;
   Mu(IGLO,IGLO)+=MuT;
   % Matrice de couplage
   % -\int\div\vvec q ou \int\vvec\cdot\grad q
   [BNCx,BNCy,BCRx,BCRy,BFSx,BFSy]=BxyNCTriangle(aire,EdgNormT,ordre);
   PGLO=NumPloc(t,:);
   Bx(PGLO,IGLO)+=BNCx;
   By(PGLO,IGLO)+=BNCy;
   Mp(PGLO,PGLO)=aire*Mploc;
   invMp(PGLO,PGLO)=(1/aire)*invMploc;
end