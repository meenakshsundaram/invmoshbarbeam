%
%    This file is a part of the invmoshbarbeam a matlab library. This 
%    is free software: you can redistribute it and/or modify it under 
%    the terms of the GNU Lesser General Public License as published by the 
%    Free Software Foundation, either version 3 of the License, or 
%    (at your option) any later version.
%
%    invmoshbarbeam is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with invmoshbarbeam.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2012 Meenakshi Sundaram   
%
%
function [U, Lambda]=fem(E,Rho,Radius,Ledis,numele,eignum,fixeddofs,tol)
%
% Function computes the linear finite element mode shape for a beam
% Input:
%   E           -   Youngs modulus of each element      (1,numele)
%   Rho         -   Density of each element             (1,numele)
%   Radius      -   Radius of each element              (1,numele)
%   Ledis       -   Length of the each element          (1,numele)
%   numele      -   Number of elements in the beam
%   eignum      -   The number of eigenmode shape desired
%   fixeddofs   -   The fixed degrees of freedom for the beam
%   tol         -   Tolerance to which the eigenvalue is determined
% Output:
%   U           -   Eigen Mode Shape normalized by maximum value
%   Lambda      -   Eigenvalue
%

    % Total degrees of freedom
    totdofs=2*(numele+1);

    % Assembly setup
    KX=zeros(16*numele,1);KY=zeros(16*numele,1);KZ=zeros(16*numele,1);
    MX=zeros(16*numele,1);MY=zeros(16*numele,1);MZ=zeros(16*numele,1);

    % Ntriplets
    ntriplets=0;
    
    % Moment of Inertia
    MI=pi*Radius.^4/4;
    
    % Area
    A=pi*Radius.^2;

    % Assembly of elements
    for ele=1:numele
    
        eledofs=[2*ele-1 2*ele 2*ele+1 2*(ele+1)];
        gdofs=[ele ele+1];
    
        [Ke1,Ke2,Me1,Me2]=elestiff(Ledis(ele));
        Keloc=Ke1*MI(gdofs(1))+Ke2*MI(gdofs(2));
        Keloc=Keloc*E(ele);
    
        Meloc=Me1*A(gdofs(1))+Me2*A(gdofs(2));
        Meloc=Meloc*Rho(ele);
    
        for i = 1:4
            for j=1:4
    
                ntriplets=ntriplets+1;
                KX(ntriplets)=eledofs(i); 
                KY(ntriplets)=eledofs(j); 
                KZ(ntriplets)=Keloc(i,j);
                MX(ntriplets)=eledofs(i); 
                MY(ntriplets)=eledofs(j); 
                MZ(ntriplets)=Meloc(i,j);
            
            end
        end
    
    end

    K=sparse(KX,KY,KZ,totdofs,totdofs);
    M=sparse(MX,MY,MZ,totdofs,totdofs);

    alldofs=1:totdofs;
    freedofs=setdiff(alldofs,fixeddofs);

    %solving
    options.tol=tol;
    options.issym=1;
    options.isreal=1;
    options.disp=0;
    U=zeros(totdofs,1);
    
    [eigvec,eigval]=eigs(K(freedofs,freedofs),M(freedofs,freedofs),eignum,0,options);

    U(freedofs)=eigvec(:,1);
    Lambda=eigval(1,1);

    %normalisation
    [~,ind]=max(abs(U));
    U=U/U(ind);