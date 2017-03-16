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
function [A,Lambda]=poweriter(psK,psM,K,M,Pow,Comul,Vol,Tol)
%
% Function computes the real eigenvalue, eigenvector using modified power
% iteration for a generalized rectangular eigen system 
%
% Input:
%   psK         -   K'*K;
%   psM         -   M'*M;
%   K,M         -   The generalized rectangular eigen system K A=lambda M A.^(1/Pow)
%   Pow         -   The power to which the multiplying vector is raised to
%   Comul       -   The vector which you can take the inner product of the
%                   area with to find total volume
%   Vol         -   The total volume of the beam
%   Tol         -   Tolerance to which the eigenvalue is determined
% Output:
%   A           -   shape parameter of each element
%   Lambda      -   eigenvalue corressponding 
%

    % Num elements
    num=size(psK,1);
    
    % Pseudo inverse
    Mat=psK\psM;

    % Initialization
    A=ones(num,1);

    % Residual and Change
    Res=1.+Tol;
    Change=1.0+Tol;

    % Power iterations
    while(Res>Tol || Change>Tol)
    
        % Storing the old shape area
        oldA=A;

        % Power iteration
        A=Mat*A;

        % Reducing the Shape parameter to the required vector (useful in beams)
        A=(abs(A).^(1/Pow));
        
        % Normalization
        A=A/max(A);
    
        % Eigenvalue
        Lambda=(A'*(A.^Pow))/(A'*Mat*A);

        % Residual
        Res=norm(K*(A.^Pow)/Lambda-M*A);
        
        % Change from old value
        Change=norm(oldA-A);    
       
    end
    
    % Rescaling to match volume
    obVol=Comul'*A;
    rati=Vol/obVol;
    A=rati*A;
    
    % Readjusting the eigenvalue
    Lambda=Lambda*(rati^(Pow-1));

