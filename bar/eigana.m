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
function [A,Lambda]=eigana(psK,psM,K,M,CoMul,Vol,tol)
%
% Function computes the real eigenvalue, eigenvector using eigenvalue
% routine - possible only when the vector involved is linear
%
% Input:
%   psK         -   K'*K;
%   psM         -   M'*M;
%   K,M         -   The generalized rectangular eigen system K A=lambda M A
%   Pow         -   The power to which the multiplying vector is raised to
%   Comul       -   The vector which you can take the inner product of the
%                   area with to find total volume
%   Vol         -   The total volume of the bar
%   Tol         -   Tolerance to which the eigenvalue is determined
% Output:
%   A           -   Area of each element
%   Lambda      -   eigenvalue corressponding 
%


    %Eigenvalue Analysis
    [eigvec,eigval]=eig(full(psK),full(psM),'chol');
    numele=size(psK,1);
    %PostProcessing
    %Finding Valid Areas
    Num=0;
    for ele=1:numele
        %the values should be reald
        if(isreal(eigvec(:,ele)) && isreal(eigval(ele,ele)))
            [val,loc]=max(abs(eigvec(:,ele)));
            %Normalization
            eigvec(:,ele)=eigvec(:,ele)/eigvec(loc,ele);
            flag=0;
            for j=1:numele
                %after normalization none of the quantities should be
                %negative
                if(eigvec(j,ele)<tol)
                    flag=-1;
                    break;
                end
            end
            if(flag==0)
                err=norm(K*eigvec(:,ele)/eigval(ele,ele)-M*eigvec(:,ele),2);
                %check the equation solving
                if(err<tol)
                    Num=Num+1;
                    A(:,Num)=eigvec(:,ele);
                    obVol=CoMul'*A;
                    A=(Vol/obVol)*A;
                    Lambda(Num)=eigval(ele,ele);
                end
            end            
        end
    end
    if(Num==0)
        A=[];
        Lambda=[];
    end