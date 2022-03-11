%{
/****************************************************************************
* Copyright (c) 2022, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% MassLGTriangle.m:
% Matrices de masse EF de Lagrange P1 et P2 (aire=1)
%
% SYNOPSIS [Mass1,invMass1,Mass01,Mass2,invMass2,Mass02,Mass12]= MassLGTriangle
%          
% OUTPUT - Mass1 : matrice de masse P1 locale
%        - invMass1 : inverse de la matrice de masse P1 locale
%        - Mass01 : matrices de couplage P0-P1
%        - Mass2  : matrice de masse P2 locale
%        - invMass2 : inverse de la matrice de masse P2 locale
%        - Mass01 : matrices de couplage P0-P2
%        - Mass12 : matrices de couplage P1-P2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mass1,invMass1,Mass01,Mass2,invMass2,Mass02,Mass12]=MassLGTriangle()
%
% MATRICES DE MASSE
%
% P1
Mass1 = (1/12)*[2 1 1;1 2 1;1 1 2];
%
invMass1 = [9 -3 -3;-3 9 -3;-3 -3 9];
%
% P2
Mass2 = (1/180)*[ 6 -1 -1 -4  0  0; ...
                 -1  6 -1  0 -4  0; ...
                 -1 -1  6  0  0 -4; ...
                 -4  0  0 32 16 16; ...
                  0 -4  0 16 32 16; ...
                  0  0 -4 16 16 32];
%
invMass2 = (1/8)*[ 288    48    48    48   -12   -12; ...
                    48   288    48   -12    48   -12; ...
                    48    48   288   -12   -12    48; ...
                    48   -12   -12    78   -27   -27; ...
                   -12    48   -12   -27    78   -27; ...
                   -12   -12    48   -27   -27    78];
%
% Couplage P0-P1
Mass01=(1/3)*[1 1 1];
% Couplage P0-P2
Mass02=(1/3)*[ -0.5 -0.5 -0.5 1 1 1];
% Couplage P1-P2
Mass12=(1/60)*[ 2,-1,-1,4,8,8;
               -1, 2,-1,8,4,8;
               -1,-1, 2,8,8,4];
