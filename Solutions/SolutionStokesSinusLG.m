%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% SolutionStokesSinusLG.m:
% Solution du probleme de Stokes "Sinus"

%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= npi*[ sin(npi*y)*{(b-2*nu*npi)cos(npi*x)+nu*npi} ; 
%                                sin(npi*x)*{(b+2*nu*npi)cos(npi*y)-nu*npi)}]
%    div U = 0;
%
% U(x,y)=[(1-cos(npi*x))*sin(npi*y) ; 
%         (cos(npi*y)-1)*sin(npi*x)]; 
% p(x,y)=sin(npi*x)*sin(npi*y);
%
%
% SYNOPSIS [Uxex,Uyex,p]=SolutionStokesSinusLG
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
% OUTPUT - [Uxex,Uyex,Pex] : la solution exacte au points de discretisation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uxex,Uyex,Pex]=SolutionStokesSinusLG(kLG)
global CoorNeu CoorMil npi
%
X=npi*CoorNeu(:,1); Y=npi*CoorNeu(:,2);
if (kLG==2)
  X=[X;npi*CoorMil(:,1)]; Y=[Y;npi*CoorMil(:,2)];
end
cosX=cos(X); cosY=cos(Y); sinX=sin(X); sinY=sin(Y);
Uxex=(1-cosX).*sinY;
Uyex=(cosY-1).*sinX;
Pex=sinX.*sinY;
%
