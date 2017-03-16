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
function [U,Lambda]=fem(E,Rho,A,L,numele,eignum,fixeddofs,tol)
%
% Function computes the linear finite element mode shape for a bar
% Input:
%   E           -   Youngs modulus of each element (1,numele)
%   Rho         -   Density of each element        (1,numele)
%   A           -   Area of each element           (1,numele)
%   L           -   Total length of the bar        
%   numele      -   Number of elements in the bar
%   eignum      -   The number of eigenmode shape desired
%   fixeddofs   -   The fixed degrees of freedom for the bar
%   tol         -   Tolerance to which the eigenvalue is determined
% Output:
%   U           -   Eigen Mode Shape normalized by maximum value
%   Lambda      -   Eigenvalue
%

    % Length of each element
    Le=L/numele;

    % Total degrees of freedom
    totdofs=2*numele+1;

    %Mass Matrices associated with each node
    Me1=[7/60 1/15 -1/60;
        1/15  4/15  0;
        -1/60 0 1/60];
    Me1=Me1*Le;

    Me2=[1/60 0 -1/60;
        0 4/15 1/15;
        -1/60 1/15 7/60];
    Me2=Me2*Le;

    %Stiffness Matrices associated with each node
    Ke1=[11/6 -2 1/6;
        -2  8/3 -2/3
        1/6 -2/3 1/2];
    Ke1=Ke1/Le;
    
    Ke2=[1/2 -2/3 1/6;
        -2/3 8/3 -2;
        1/6 -2 11/6];
    Ke2=Ke2/Le;

    %Assembly
    KX=zeros(9*numele,1);KY=zeros(9*numele,1);KZ=zeros(9*numele,1);
    MX=zeros(9*numele,1);MY=zeros(9*numele,1);MZ=zeros(9*numele,1);
    
    ntriplets=0;
    for ele=1:numele    
    
        eledofs=[2*ele-1 2*ele 2*ele+1];
        Adofs=[ele ele+1];
    
        Ke=Ke1*A(Adofs(1))+Ke2*A(Adofs(2));
    
        Me=Me1*A(Adofs(1))+Me2*A(Adofs(2));
    
        Ke=Ke*E(ele);
    
        Me=Me*Rho(ele);
    
        for i = 1:3
    
            for j=1:3
        
                ntriplets=ntriplets+1;
                KX(ntriplets)=eledofs(i); KY(ntriplets)=eledofs(j); KZ(ntriplets)=Ke(i,j);
                MX(ntriplets)=eledofs(i); MY(ntriplets)=eledofs(j); MZ(ntriplets)=Me(i,j);
            
            end
        end
    
    end

    K=sparse(KX,KY,KZ,totdofs,totdofs);
    M=sparse(MX,MY,MZ,totdofs,totdofs);

    % Boundary Conditions
    alldofs=1:totdofs;
    freedofs=setdiff(alldofs,fixeddofs);

    %solving for the eigenvalue, eigenmodeshape
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