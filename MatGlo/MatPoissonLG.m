%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% MatPoissonLG.m:
% Assemblages des matrices de masse et de raideur
%
% avec les elements finis de Galerkin continus (polynomes de Lagrange P1 ou P2)
%
% SYNOPSIS [Ku,Mu] = MatPoissonLG(CoorNeu,NumTri,TriEdg,EdgTri,Lga2,EdgNorm,Aires,ordre,nuT)
%          
% INPUT  - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(NbEdg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - Lga2(NbEdg,1)    : longueurs des aretes au carre
%        - EdgNorm(NbEdg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1)   : aires des triangles
%        - ordre : 1 ou 2
%        - nuT : coefficient de raideur
%
% OUTPUT - Ku : matrice de raideur
%        - Mu : matrice de masse 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ku,Mu] = MatPoissonLG(CoorNeu,NumTri,TriEdg,EdgTri,Lga2,EdgNorm,Aires,ordre,nuT)
%
Nbpt=size(CoorNeu,1);
Nbtri=size(NumTri,1);
%
% Matrices de masse locales pour aire=1 !!
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle;

if (ordre==1)
  Ndof_u=Nbpt;
  Numloc=NumTri; 
  MLG=MT1;
end
%
if (ordre==2)
  NbEdg=size(EdgTri,1);
  Ndof_u=Nbpt+NbEdg;
  Numloc=[NumTri,TriEdg+Nbpt];
  MLG=MT2;
end
Ku=sparse(Ndof_u,Ndof_u); Mu=Ku;
for t=1:Nbtri
  aire=Aires(t);
  AGLO=TriEdg(t,:);
  EdgNormT=EdgNorm(AGLO,:);
  Lga2T=Lga2(AGLO);
  for iloc=1:2
    if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
    end
  end
  EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
  %Lga2T(3)=dot(EdgNormT(3,:),EdgNormT(3,:));
  % On suppose nu constant par element
  [KLG,K1,K2]=KriLGTriangle(aire,EdgNormT,Lga2T,ordre);
  KuT=nuT(t)*KLG; MuT=aire*MLG;
  %% Ku et Mu
  IGLO=Numloc(t,:);
  Ku(IGLO,IGLO)+=KuT;
  Mu(IGLO,IGLO)+=MuT;
end % end_for t=1:Nbtri
%

