%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% MatMassLG.m:
% Assemblage matrice de masse EF de Lagrange P1 ou P2
%
% SYNOPSIS Mu = MatMassLG
%          
% GLOBAL - Nbpt : nombre de sommets
%        - Nbedg : nombre d'aretes
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'Edg opposee au sommet NumTri(l,i)
%                  (3 numeros des Edg - matrice entiere Nbtri x 3)
%		     - Aires(Nbtri,1) : aires des triangles
% INPUT  - k : 0, 1 ou 2
%
% OUTPUT - Mu : matrice de masse 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mu = MatMassLG(kLG)
%
global Nbpt
global Nbtri Aires NumTri
global Nbedg TriEdg
%
% Matrices de masse locales pour aire=1 !!
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle();
%
% ORDRE 0
if (kLG==0)
  Mu=diag(sparse(Aires),0);
else
  % ORDRE 1
  if (kLG==1)
     MT=MT1; Ndof=Nbpt; nloc=3; Numloc=NumTri;
  end
  % ORDRE 2
  if (kLG==2)
     MT=MT2; Ndof=Nbedg+Nbpt; nloc=6; Numloc=[NumTri,Nbpt+TriEdg];
  end
  Mu=sparse(Ndof,Ndof);
  for t=1:Nbtri
     Mu_loc=Aires(t)*MT;
     IGLO=Numloc(t,:);
     Mu(IGLO,IGLO)+=Mu_loc;
  end % end_for t=1:Nbtri
end
