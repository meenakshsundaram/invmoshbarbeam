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
function [Ke1,Ke2,Me1,Me2]=elestiff(Le)
%
% Function computes the element stiffness
% Input:
%   Le          -   Length of the element
% Output:
%   Ke1,Ke2     -   Element stiffness matrices   
%   Me1,Me2     -   Element mass matrices
%

%Mass Matrix

    Me1	=[	(2*Le)/7	Le^2/28		(9*Le)/140	-(Le^2/60);
        	Le^2/28		Le^3/168	Le^2/70		-(Le^3/280);
            (9*Le)/140	Le^2/70		(3*Le)/35	-(Le^2/60);
            -(Le^2/60)	-(Le^3/280)	-(Le^2/60)	Le^3/280;];

    Me2	=[	(3*Le)/35	Le^2/60		(9*Le)/140	-(Le^2/70);
        	Le^2/60		Le^3/280	Le^2/60		-(Le^3/280);
            (9*Le)/140	Le^2/60		(2*Le)/7	-(Le^2/28);
            -(Le^2/70)	-(Le^3/280)	-(Le^2/28)	Le^3/168;];

%Stiffness Matrix

    Ke1	=[	6/Le^3 		4/Le^2 		-(6/Le^3) 	2/Le^2;
        	4/Le^2		3/Le		-(4/Le^2)	1/Le;
            -(6/Le^3)	-(4/Le^2)	6/Le^3		-(2/Le^2);
            2/Le^2		1/Le		-(2/Le^2)	1/Le;];

    Ke2	=[	6/Le^3		2/Le^2		-(6/Le^3)	4/Le^2;
        	2/Le^2		1/Le		-(2/Le^2)	1/Le;
            -(6/Le^3)	-(2/Le^2)	6/Le^3		-(4/Le^2);
            4/Le^2		1/Le		-(4/Le^2)	3/Le;];