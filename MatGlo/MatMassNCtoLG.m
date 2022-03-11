%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatMassNCtoLG.m:
% Matrice de projection NC vers LG P1 ou P2
%
% SYNOPSIS MatP=MatMassNCtoLG(Nbpt,NbEdg,NumTri,TriEdg,Aires,ordre)
%          
% INPUT  - Nbpt : nombre de sommets
%        - NbEdg : nombre d'aretes
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'Edg opposee au sommet NumTri(l,i)
%                  (3 numeros des Edg - matrice entiere Nbtri x 3)
%		     - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2 
%
% OUTPUT - MatP : Matrice de projection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatP=MatMassNCtoLG(Nbpt,NbEdg,NumTri,TriEdg,Aires,ordre)

Nbtri=size(NumTri,1);
%
if (ordre==1)
  MatP=sparse(Nbpt,NbEdg);
  for t=1:Nbtri
    aire6=Aires(t)/6;
    for iloc=1:3
      iglo=NumTri(t,iloc);
      for aloc=1:3
        if (aloc~=iloc)
           aglo=TriEdg(t,aloc);
           MatP(iglo,aglo)+=aire6;
        endif       
      endfor
    endfor
  end
endif

if (ordre==2)
  Ndof=Nbpt+NbEdg;
  MatP=[Mu_LG,sparse(Ndof,Nbtri)];
  TGLO=Ndof;
  for t=1:Nbtri
    aire5=Aires(t)/5;
    aire30=Aires(t)/30;
    TGLO=TGLO+1;
    IGLO=NumTri(t,:);
    MatP(IGLO,TGLO)=aire30*ones(3,1);
    AGLO=TriEdg(t,:)+Nbpt;
    MatP(AGLO,TGLO)=aire5*ones(3,1);
  end
endif
