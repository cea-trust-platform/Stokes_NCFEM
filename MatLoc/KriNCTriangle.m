%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% KriNCTriangle.m:
% Matrices de raideur interne EF de CR P1-NC (ordre=1) ou FS (ordre=2) 
%                             pour un triangle donne
% 
% Pour tester : aire=0.5; EdgNorm=[1,1;-1,0;0,-1]; Lga2=[2,1,1];
%
% SYNOPSIS [KNC,KCR,KFS] = KriNCTriangle(aire,EdgNorm,Lga2,ordre)
%          
% INPUT  - aire = aire du triangle
%        - EdgNorm(3,2) : coordonnees (x, y) des faces normales
%        - Lga2(3,1)    : longueurs des aretes au carre
%        - ordre            : ordre d'approximation
% OUTPUT - KNC : matrice de raideur locales pour EF de NC de l'ordre concerne
%        - KCR : matrice de raideur locales pour EF de CR P1-NC
%        - KFS : matrice de raideur locales pour EF de FS P2+BulleTriangle 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KNC,KCR,KFS]=KriNCTriangle(aire,EdgNorm,Lga2,ordre)
%
% MATRICES DE RAIDEUR
%
[KLG,K1,K2]=KriLGTriangle(aire,EdgNorm,Lga2,ordre);
KCR=4*K1;
KFS=[];
if (ordre==1)
  KNC=KCR;
end
if (ordre==2)
   Ktt=3*sum(Lga2)/(4*aire);
   K2t=zeros(6,1);
   for i=1:3        
       j=mod(i,3)+1;
       k=6-(i+j);
       K2t(i)=-0.5*KCR(i,i);
       K2t(i+3)=-KCR(j,k);
   endfor   
   KFS=[K2,K2t;K2t',Ktt];
   KNC=KFS;
   % CHECK
##   XY=[0,0;1,0;0,1];
##   [xyp,wp,lambda,np]=IntTri_Ham7(XY);
##   Glambda=-(0.5/aire)*EdgNorm;
##   wp=wp*aire;
##   GphiT=-6*lambda*Glambda;
##   ij=[2,3,1];
##   for i=1:3
##     GphiS=(4*lambda(:,i)-1)*Glambda(i,:);
##     j=ij(i); k=6-(i+j);
##     K2t(i)=wp*(sum((GphiT.*GphiS),2));
##     GphiM=4*(lambda(:,j)*Glambda(k,:)+lambda(:,k)*Glambda(j,:));
##     K2t(i+3)=wp*(sum((GphiT.*GphiM),2));
##   endfor
##   Ktt=wp*(sum((GphiT.*GphiT),2));
##   KFS2=[K2,K2t;K2t',Ktt]
endif