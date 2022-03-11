%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% KriLGTriangle.m:
% Matrices de raideur interne EF de Lagrange P1 et P2 pour un triangle donne
% 
% Pour tester : aire=0.5; EdgNorm=[1,1;-1,0;0,-1]; Lga2=[2,1,1]
%
% SYNOPSIS [KLG,K1,K2] = KriLGTriangle(aire,EdgNorm,Lga2,ordre)
%          
% INPUT  - aire : aire du triangle
%        - EdgNorm(3,2) : coordonnees (x, y) des faces normales
%        - Lga2(3) : longueurs des aretes au carre
%        - ordre            : ordre d'approximation
% OUTPUT - KLG, K1, K2 : matrices de raideur locales pour EF de Lagrange ordre, P1 P2
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KLG,K1,K2]=KriLGTriangle(aire,EdgNorm,Lga2,ordre)
%
% MATRICES DE RAIDEUR
%
% P1
%
% K1 = Produit Scalaire des gradients des coordonnees barycentriques
%      
%
K1=zeros(3,3);
a4=1/(4*aire);
for i=1:3
  K1(i,i)=a4*Lga2(i);
end
K1(1,2)=a4*dot(EdgNorm(1,:),EdgNorm(2,:));
K1(2,1)=K1(1,2);
K1(1,3)=-K1(1,1)-K1(1,2);
K1(3,1)=K1(1,3);
K1(2,3)=-K1(2,2)-K1(2,1);
K1(3,2)=K1(2,3);

if (ordre==1)
   KLG=K1;
   K2=[];
end
if (ordre==2)
   K2=zeros(6,6);
   for i=1:3
     K2(i,i)=K1(i,i);
     for j=i+1:3
       K2(i,j)=-K1(i,j)/3; K2(j,i)=K2(i,j);
     end
     j=mod(i,3)+1;
     k=6-(i+j);       
     K2(i+3,i+3)=8*(K1(i,i)-K1(j,k))/3;
   end
   K2(1,5)=4*K1(1,3)/3; K2(5,1)=K2(1,5);
   K2(1,6)=4*K1(1,2)/3; K2(6,1)=K2(1,6);
   K2(2,4)=4*K1(2,3)/3; K2(4,2)=K2(2,4);
   K2(2,6)=K2(1,6);     K2(6,2)=K2(2,6);
   K2(3,4)=K2(2,4);     K2(4,3)=K2(3,4);
   K2(3,5)=K2(1,5);     K2(5,3)=K2(3,5);
   K2(4,5)=8*K1(1,2)/3; K2(5,4)=K2(4,5);
   K2(4,6)=2*K2(3,5);   K2(6,4)=K2(4,6);
   K2(5,6)=2*K2(3,4);   K2(6,5)=K2(5,6);
##   verif=0;
##   if (verif==1)
##      Un = ones(6,1);
##      KLGn=K2*Un; printf(' ordre 2, max Kri*Un = %7.2e\n',max(abs(KLGn)));
##   end
   KLG=K2;
end