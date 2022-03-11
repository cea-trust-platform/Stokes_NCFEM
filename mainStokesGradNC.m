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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% mainStokesGradNC.m:
%
% Convergence en maillage, EF non conformes CR ou FS 
%   Pb de Stokes dans un carre [0,1]*[0,1] :
%   -nu*Delta U + grad p= 3[x^2,y^2]
%    div U = 0;
%
% U(x,y)=0
% p(x,y)=x^3+y^3-1/2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
m0=2; ef0=2; 
m1=2; ef1=2;
%
global nu=1;
global RT=0;
global algOct=1;
%
nmesh=m1-m0+1; nef=ef1-ef0+1;
meshstep=[0.1,0.05,0.025,0.0125,0.00625];
TabNbpt  =zeros(1,nmesh);%[142, 568,2212, 8558, 34239];
TabNbtri =zeros(1,nmesh);%[242,1054,4262,16794, 68476];
TabNbedg =zeros(1,nmesh);%[383,1621,6473,25351,102074];
TabNFS   =zeros(1,nmesh);%[383,1621,6473,25351,102074];
Eu0=zeros(nef,nmesh); Eu1=zeros(nef,nmesh); EuDiv=zeros(nef,nmesh); Ep0=zeros(nef,nmesh);
im=0; 
%
% Pas d'image, que le graphe
%fig=-nmesh*nef;
%
fig=-nmesh*nef+1;
for mii=m0:m1
  im=im+1;
  fprintf(' Maillage Square_h%i\n',mii);
  filename =sprintf('Square_h%i',mii);
  global Nbpt CoorNeu CoorNeu2 RefNeu
  global Nbtri CoorBary Aires NumTri NumTri2 TriEdg
  global Nbedg CoorMil RefEdg Lga2 EdgNorm EdgTri
  global ordre mi num
  mi=mii; num=30;
  %
  [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,NumEdgB,CoorMil,...
         RefEdg,RefEdgB,TriEdg,EdgTri,SomOpp,Lga2,EdgNorm,Aires]=readmeshfiles(filename);
  Nbpt = size(CoorNeu,1); TabNbpt(im) = Nbpt;
  Nbtri= size(NumTri,1) ; TabNbtri(im)= Nbtri;
  Nbedg= size(NumEdg,1) ; TabNbedg(im)= Nbedg; 
  TabNFS(im)=Nbpt+Nbedg+Nbtri;
  %
  efi=0;
  for ef=ef0:ef1
    fig=fig+1;
    efi=efi+1;
    ordre=ef;
    fprintf('Ordre P%i, RT=%i, nu=%7.2e.\n',ordre,RT,nu);
    [Eu0(efi,im),Eu1(efi,im),Ep0(efi,im),fig]=SolveStokesGradNC(fig);
  end
end
%
if (nmesh>1)
  logH =log10(meshstep);   dlogH =logH (1:nmesh-1)-logH (2:nmesh);
  logN1=0.5*log10(TabNbpt);   dlogN1=logN1(1:nmesh-1)-logN1(2:nmesh);
  % vitesse
  logU0=log10(Eu0); dlogU0=logU0(:,1:nmesh-1)-logU0(:,2:nmesh);
  logU1=log10(Eu1); dlogU1=logU1(:,1:nmesh-1)-logU1(:,2:nmesh);
  % Pression
  logP0=log10(Ep0); dlogP0=logP0(:,1:nmesh-1)-logP0(:,2:nmesh);
  %
  invdlogH=dlogH.^(-1);
  invdlogN=dlogN1.^(-1);
  %
  % Convergence en pas du maillage
  vtauU0=zeros(nef,nmesh-1); tauU0=zeros(nef,1);
  vtauU1=zeros(nef,nmesh-1); tauU1=zeros(nef,1);  
  vtauP0=zeros(nef,nmesh-1); tauP0=zeros(nef,1);
  % Convergence en nombre de ddl
  vtauU0P=zeros(nef,nmesh-1); tauU0P=zeros(nef,1); 
  vtauU1P=zeros(nef,nmesh-1); tauU1P=zeros(nef,1); 
  vtauP0P=zeros(nef,nmesh-1); tauP0P=zeros(nef,1); 
  %
  for i=1:nef
    % Convergence en pas du maillage
    vtauU0(i,:)=dlogU0(i,:).*invdlogH; tauU0(i)=sum(vtauU0(i,:))/(nmesh-1);
    vtauU1(i,:)=dlogU1(i,:).*invdlogH; tauU1(i)=sum(vtauU1(i,:))/(nmesh-1);
    vtauP0(i,:)=dlogP0(i,:).*invdlogH; tauP0(i)=sum(vtauP0(i,:))/(nmesh-1);
    % Convergence en nombre de ddl
    vtauU0P(i,:)=-dlogU0(i,:).*invdlogN; tauU0P(i)=sum(vtauU0P(i,:))/(nmesh-1);
    vtauU1P(i,:)=-dlogU1(i,:).*invdlogN; tauU1P(i)=sum(vtauU1P(i,:))/(nmesh-1);
    vtauP0P(i,:)=-dlogP0(i,:).*invdlogN; tauP0P(i)=sum(vtauP0P(i,:))/(nmesh-1); 
  end
  if (fig>=0)
    efi=0;
    for ef=ef0:ef1
      efi=efi+1;
      titre=sprintf('Convergence EF NC, ef %i',ef);      
      fig=fig+1;
      figure(fig);
      plot(logH(1,m0:m1),logU0(efi,:),'-b;EU_0;',...
           logH(1,m0:m1),logU1(efi,:),'-r;EU_1;',...
           logH(1,m0:m1),logP0(efi,:),'-g;EP_0;');
      xlabel ('log10(h)');
      ylabel ('log10(Erreur)');
      title (titre);
    endfor
  endif
endif
%
if (nef==2)
   entete=sprintf('eU_0 CR  eU_1 CR  eP_0 ef%i  eU_0 FS  eU_1 FS  eP_0 ef%i  h      NCR\n',ef0,ef1);
   TabNinc=[TabNbedg',TabNFS'];
   if (RT==0)
     sortie=sprintf('StokesGrad-nu=%7.2e',nu);
   else
     sortie=sprintf('StokesGrad-RT-nu=%7.2e',nu);
   end
end
if (nef==1)
    if (ef0==1)
      entete=sprintf('eU_0 CR  eU_1 CR  eP_0 ef%i  h          NCR\n',ef0);
      if (RT==0)
       sortie=sprintf('StokesGradCR-nu=%7.2e',nu);
      else
       sortie=sprintf('StokesGradCR-RT-nu=%7.2e',nu);
      end
      TabNinc=TabNbedg';
    else
      entete=sprintf('eU_0 FS  eU_1 FS  eP_0 ef%i  h          NFS\n',ef0);
      if (RT==0)
       sortie=sprintf('StokesGradFS-nu=%7.2e',nu);
      else
       sortie=sprintf('StokesGradFS-RT-nu=%7.2e',nu);
      end
      TabNinc=TabNFS';
    end     
end
fid=fopen(sortie,'w');    
fprintf(fid,entete);
for j=1:nmesh
  for i=1:nef
      fprintf(fid,'%7.2e %7.2e %7.2e ',Eu0(i,j),Eu1(i,j),Ep0(i,j));
  end
  fprintf(fid,'%7.2e %5i %5i',meshstep(1,j),TabNinc(j,:));
  fprintf(fid,'% i\n',j);
end
if (nmesh>1)
  for i=1:nef
      fprintf(fid,'%7.2e %7.2e  %7.2e ',tauU0(i),tauU1(i),tauP0(i));
  end
  fprintf(fid,'\n');
  for i=1:nef
      fprintf(fid,'%7.2e %7.2e  %7.2e ',tauU0P(i),tauU1P(i),tauP0P(i));
  end
  fprintf(fid,'\n'); 
end
fclose(fid);
