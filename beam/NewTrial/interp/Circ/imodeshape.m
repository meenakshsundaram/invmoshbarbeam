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
function [Radius,Lambda]=imodeshape(E,Rho,U,Ledis,numele,fixeddofs,Vol,tol)
%
% Function computes the linear finite element radius
% for a given mode shape of the beam
% Input:
%   E           -   Youngs modulus of each element       (1,numele)
%   Rho         -   Density of each element              (1,numele)
%   U           -   The Displacement,slope of each node  (1,2*(numele+1))
%   Ledis       -   Length of each element               (1,numele)
%   numele      -   Number of elements in the beam
%   fixeddofs   -   The fixed degrees of freedom for the beam
%   Vol         -   The total volume of the beam
%   tol         -   Tolerance to which the eigenvalue is determined
% Output:
%   Radius      -   Radius of each element
%   Lambda      -   eigenvalue corressponding to the modeshape
%
    % Total degrees of freedom
    totdofs=2*(numele+1);

    %Assembly
    KX=zeros(8*numele,1);KY=zeros(8*numele,1);KZ=zeros(8*numele,1);
    MX=zeros(8*numele,1);MY=zeros(8*numele,1);MZ=zeros(8*numele,1);

    % Assembling elements
    ntriplets=0;
    for ele=1:numele    
    
        eledofs=2*ele-1:2*(ele+1);
        gdofs=[ele ele+1];
        [Ke1,Ke2,Me1,Me2]=elestiff(Ledis(ele));
        
        Kvec1=E(ele)*Ke1*U(eledofs)*pi/4;
        Kvec2=E(ele)*Ke2*U(eledofs)*pi/4;
        
        Mvec1=Rho(ele)*Me1*U(eledofs)*pi;
        Mvec2=Rho(ele)*Me2*U(eledofs)*pi;
    
        for i = 1:4
            	
            ntriplets=ntriplets+1;
                
            KX(ntriplets)=eledofs(i);
            KY(ntriplets)=gdofs(1);
            KZ(ntriplets)=Kvec1(i);
            	
            MX(ntriplets)=eledofs(i);
            MY(ntriplets)=gdofs(1);
            MZ(ntriplets)=Mvec1(i);
            
        end
        
        for i = 1:4
            
            ntriplets=ntriplets+1;
            
            KX(ntriplets)=eledofs(i);
            KY(ntriplets)=gdofs(2);
            KZ(ntriplets)=Kvec2(i);
            
            MX(ntriplets)=eledofs(i);
            MY(ntriplets)=gdofs(2);
            MZ(ntriplets)=Mvec2(i);
            
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

    CoMul=[Ledis(1)/2;(Ledis(1:end-1)+Ledis(2:end))/2;Ledis(end)/2]*pi;
    [sqRadius,Lambda]=poweriter(psK,psM,K(freedofs,:),M(freedofs,:),2,CoMul,Vol,tol);

    Radius=sqrt(sqRadius);
