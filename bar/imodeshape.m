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
function [A,Lambda]=imodeshape(E,Rho,U,L,numele,fixeddofs,Vol,tol,flag)
%
% Function computes the linear finite element area 
% for a given mode shape of the bar
% Input:
%   E           -   Youngs modulus of each element (1,numele)
%   Rho         -   Density of each element        (1,numele)
%   U           -   The Displacement of each node  (1,numele+1)
%   L           -   Total length of the bar        
%   numele      -   Number of elements in the bar
%   fixeddofs   -   The fixed degrees of freedom for the bar
%   Vol         -   The total volume of the bar
%   tol         -   Tolerance to which the eigenvalue is determined
%   flag        -   The flag which dictates which solver will be used 
%                   the iterative (1) or the eigen solver (anything else)
% Output:
%   A           -   Area of each element
%   Lambda      -   eigenvalue corressponding 
%


    Le=L/numele;
    totdofs=2*numele+1;

    %Mass Matrix
    Me1=[7/60 1/15 -1/60;
        1/15  4/15  0;
        -1/60 0 1/60];
    Me1=Me1*Le;

    Me2=[1/60 0 -1/60;
        0 4/15 1/15;
        -1/60 1/15 7/60];
    Me2=Me2*Le;

    %Stiffness Matrix
    Ke1=[11/6 -2 1/6;
        -2  8/3 -2/3
        1/6 -2/3 1/2];
    Ke1=Ke1/Le;
    
    Ke2=[1/2 -2/3 1/6;
        -2/3 8/3 -2;
        1/6 -2 11/6];
    Ke2=Ke2/Le;

    %Assembly
    KX=zeros(6*numele,1);KY=zeros(6*numele,1);KZ=zeros(6*numele,1);
    MX=zeros(6*numele,1);MY=zeros(6*numele,1);MZ=zeros(6*numele,1);

    ntriplets=0;
    for ele=1:numele    
        eledofs=[2*ele-1 2*ele 2*ele+1];
        Adofs=[ele ele+1];    
        Kvec1=Ke1*E(ele)*U(eledofs);
        Kvec2=Ke2*E(ele)*U(eledofs);
        Mvec1=Rho(ele)*Me1*U(eledofs);
        Mvec2=Rho(ele)*Me2*U(eledofs);
        for i = 1:3
            ntriplets=ntriplets+1;
            KX(ntriplets)=eledofs(i);KY(ntriplets)=Adofs(1);KZ(ntriplets)=Kvec1(i);
            MX(ntriplets)=eledofs(i);MY(ntriplets)=Adofs(1);MZ(ntriplets)=Mvec1(i);        
        end
        for i = 1:3
            ntriplets=ntriplets+1;
            KX(ntriplets)=eledofs(i);KY(ntriplets)=Adofs(2);KZ(ntriplets)=Kvec2(i);
            MX(ntriplets)=eledofs(i);MY(ntriplets)=Adofs(2);MZ(ntriplets)=Mvec2(i);        
        end    
    end

    K=sparse(KX,KY,KZ,totdofs,numele+1);
    M=sparse(MX,MY,MZ,totdofs,numele+1);
    
    %Setting up Boundary Conditions
    alldofs=1:totdofs;
    freedofs=setdiff(alldofs,fixeddofs);

    %Pseudo Inverse
    psK=K(freedofs,:)'*K(freedofs,:);
    psM=K(freedofs,:)'*M(freedofs,:);
    CoMul=[Le/2; ones(numele-1,1)*Le; Le/2];
    if(flag==1)
        % Iterative method
        [A,Lambda]=poweriter(psK,psM,K(freedofs,:),M(freedofs,:),1,CoMul,Vol,tol);
    else
        % Eigenvalue method
        [A,Lambda]=eigana(psK,psM,K(freedofs,:),M(freedofs,:),CoMul,Vol,tol);
    end

