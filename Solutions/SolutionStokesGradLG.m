%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% SolutionStokesGradLG.m:
% Solution du probleme de Stokes "Grad"
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= 3[x^2,y^2]
%    div U = 0;
%
% U(x,y)=0
% p(x,y)=x^3+y^3-1/2;
%
% SYNOPSIS Pex=SolutionStokesGradLG
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets
% OUTPUT - Pex : la solution exacte au points de discretisation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pex=SolutionStokesGradLG( )
global CoorNeu 
%
X=CoorNeu(:,1); Y=CoorNeu(:,2);

%
