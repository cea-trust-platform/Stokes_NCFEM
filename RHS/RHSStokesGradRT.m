%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% RHSStokesGradRT.m:
%
% Calcul du second membre pour resoudre le probleme de Stokes "Grad" 
%   avec les elements finis non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= 3[x^2,y^2]
%    div U = 0;
%
% U(x,y)=0
% p(x,y)=x^3+y^3-1/2;
%
% SYNOPSIS [RHS_Ux,RHS_Uy]=RHSStokesGradRT(Mu)
%          
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
%        - RefNeu(Nbpt,1) : reference des sommets
%        - NumTri(Nbtri,3) : liste de triangles (3 numeros de sommets)
%        - Nbedg : nombre d'aretes
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - RefEdg(Nbedg,1) : Reference de chaque arete
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre : 1 ou 2
%
% OUTPUT - RHS_Ux : second membre composante x
%        - RHS_Uy : second membre composante y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHS_Ux,RHS_Uy]=RHSStokesGradRT()
%
global Nbpt CoorNeu RefNeu
global Nbtri Aires NumTri
global Nbedg TriEdg RefEdg EdgTri EdgNorm
global ordre
%
Ndof=Nbedg;
if (ordre==2)
  NP2=Nbedg+Nbpt;
  Ndof=NP2+Nbtri;
end
RHS_Ux=zeros(Ndof,1); RHS_Uy=zeros(Ndof,1);
%
% EF de CR P1-NC (ordre=1) 
if (ordre==1)
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    CoorNeuT=CoorNeu(IGLO,:);
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    awp=1.5*wp';
    FxLOC = awp.*xyp(:,1).^2; 
    FyLOC = awp.*xyp(:,2).^2;
    % (d|T|)^{-1}\int_T(x-OSi)\cdot\Fvec
    RHS_tmp=(ones(3,1)*xyp(:,1)'-CoorNeuT(:,1)*ones(1,np))*FxLOC+...
            (ones(3,1)*xyp(:,2)'-CoorNeuT(:,2)*ones(1,np))*FyLOC;
    %
    AGLO=TriEdg(t,:); 
    % vecteurs face-normale
    AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
    for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
    RHS_Ux(AGLO)+=EdgNormT(:,1).*RHS_tmp;
    RHS_Uy(AGLO)+=EdgNormT(:,2).*RHS_tmp;
    %
  end
end
%
% EF de FS (ordre=2) 
if (ordre==2)
  ILOC=[1,2,3]; JLOC=[2,3,1]; KLOC=[3,1,2];
  phiM=zeros(7,3);
  for t=1:Nbtri
    IGLO=NumTri(t,:);
    CoorNeuT =CoorNeu(IGLO,:);
    % vecteurs face-normale
    AGLO=TriEdg(t,:); EdgNormT=EdgNorm(AGLO,:);
    for iloc=1:2
     if (EdgTri(AGLO(iloc),1)~=t)
      EdgNormT(iloc,:)=-EdgNormT(iloc,:);
     end
    end
    EdgNormT(3,:)=-EdgNormT(1,:)-EdgNormT(2,:);
    % volume of the element
    aire=Aires(t);
    % Integration points
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    awp=3*aire*wp';
    FxLOC = awp.*xyp(:,1).^2; 
    FyLOC = awp.*xyp(:,2).^2;
    FLOC = [FxLOC,FyLOC];
    %
    [PiRTxx,PiRTxy,PiRTyx,PiRTyy]=PiRT1(aire,CoorNeuT,EdgNormT,lambda);
    %
    UGLO=[IGLO,AGLO+Nbpt,NP2+t];
    RHS_Ux(UGLO)+=(FxLOC'*PiRTxx+FyLOC'*PiRTxy)';
    RHS_Uy(UGLO)+=(FxLOC'*PiRTyx+FyLOC'*PiRTyy)';
    % 
  end
end
%
%% Traitement des CL de Dirichlet
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
Nbord=size(Dof_Bord,1);
RHS_Ux(Dof_Bord,1)=zeros(Nbord,1);
RHS_Uy(Dof_Bord,1)=zeros(Nbord,1);
%