%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyright : Erell Jamelot CEA
%
% BKBp.m:
%
% Calcule (Bmat*Kmat^{-1}*Bmat^T)*p avec Kmat=Rmat'*Rmat
%
% SYNOPSIS Ap=BKBtp(Rmat,Bmat,p)
%
% INPUT  - Rmat(Nu,Nu) matrice issue de la factorisation de Cholesky
%                      de la matrice de raideur de la vitesse
%        - Bmat(Np,Nu) matrice de couplage vitesse-pression
%        - Ph(Np,1)   pression
%          
% OUTPUT - Ap(np,1):=(Bmat*Kmat^{-1}*Bmat^T)*p
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ap=BKBtp(Rmat,Bmat,Ph)
  fvec=Bmat'*Ph;
  % 
  % algorithme de descente-remontée
  Uvec=Rmat\(Rmat'\fvec);
  Ap=Bmat*Uvec;
end