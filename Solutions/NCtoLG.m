%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% NCtoLG.m:
%
% Projection NC vers LG P1 ou P2
%
% SYNOPSIS [Uxh_LG,Uyh_LG,Ph_LG]=NCtoLG(Mu_LG,Uxh,Uyh,Ph)
%          
% GLOBAL - Nbpt : nombre de sommets
%        - Nbedg : nombre d'aretes
%        - Uh(Ndof,1) : vecteur a  projeter
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'Edg opposee au sommet NumTri(l,i)
%                  (3 numeros des Edg - matrice entiere Nbtri x 3)
%		     - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2 
%
% INPUT  - Uxh,Uyh,Ph
%
% OUTPUT - Uh_LG : vecteur projete
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxh_LG,Uyh_LG,Ph_LG]=NCtoLG(Uxh,Uyh,Ph)
%
global Nbpt Nbedg Nbtri NumTri TriEdg Aires ordre
%
  Mk_LG  = MatMassLG(ordre);
  ndofP=size(Ph,1);
  ordreP = (ordre-1)*ones(ndofP,1);
  Ph_LG  = DGPktoLG(Mk_LG,Ph,ordreP,ordre);
  if (ordre==1)
    MP=sparse(Nbpt,Nbedg);
    for t=1:Nbtri
      aire6=Aires(t)/6;
      for iloc=1:3
        iglo=NumTri(t,iloc);
        for aloc=1:3
          if (aloc~=iloc)
             aglo=TriEdg(t,aloc);
             MP(iglo,aglo)=MP(iglo,aglo)+aire6;
          endif       
        endfor
      endfor
    end
    PUxh=MP*Uxh; Uxh_LG=Mk_LG\PUxh;
    PUyh=MP*Uyh; Uyh_LG=Mk_LG\PUyh;
  endif

  if (ordre==2)
    [MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle;
    [MCR,MFS]=MassNCTriangle(ordre,MT2);
    Ndof=Nbpt+Nbedg;
    MP=[Mk_LG,sparse(Ndof,Nbtri)];
    TGLO=Ndof;
    for t=1:Nbtri
      aire5=Aires(t)/5;
      aire30=Aires(t)/30;
      TGLO=TGLO+1;
      IGLO=NumTri(t,:);
      MP(IGLO,TGLO)=aire30*ones(3,1);
      AGLO=TriEdg(t,:)+Nbpt;
      MP(AGLO,TGLO)=aire5*ones(3,1);
    end
    PUxh=MP*Uxh; Uxh_LG=Mk_LG\PUxh;
    PUyh=MP*Uyh; Uyh_LG=Mk_LG\PUyh;
  endif
