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
% IntTri_Ham7.m:
% Points d'integration dans le triangle de sommets de coordonnees XY(3,2)
% Formule d'integration dans un triangle avec interpolation a 7 points
%         exacte pour les polynomes d'ordre 5
% Biblio : formule d'Hammer-Marlow-Stroud, Dautray-Lions, tome 6, page 781
% INPUT  - XY(3,2) : coordonnees des sommets du triangle
% OUTPUT - xyp(7,2) : coordonnees des points d'interpolation
%        - wp(1,7) : poids associes
%        - lambda(7,3) : coordonnees barycentriques associees aux points d'interpolation
%        - np = 7 : nombre de points d'interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xyp,wp,lambda,np]=IntTri_Ham7(XY)


np=7;

c=sqrt(15);
pM=(155-c)/1200;
pP=(155+c)/1200;
wp=[9./40,pM,pM,pM,pP,pP,pP];
wp(1,1)=1-sum(wp(1,2:np));
lM=(6-c)/21; lM3=1-2*lM;
lP=(6+c)/21; lP3=1-2*lP;
lambda=[1./3,1./3,1-2/3;
        lM ,lM ,lM3;
        lM ,lM3,lM ;
        lM3,lM ,lM ;
        lP ,lP ,lP3;
        lP ,lP3,lP ;
        lP3,lP ,lP ];
xyp=lambda*XY;
