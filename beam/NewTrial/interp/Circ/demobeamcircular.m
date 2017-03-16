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
%% Initialization 
clc;
clear;

% Number of elements
numele=32;

% The Mode shape number
eignum=2;

% Fixed degrees of freedom
fixeddofs=[1,2,2*(numele)+1,2*(numele+1)];

% Total length of the beam
L=1;

% nodal points
x=linspace(0,L,numele+1)';

% Length of each element
Ledis=x(2:end)-x(1:end-1);

% Synthesized radius
Radius=sin(pi*x)*5+10;

% Area of the beam
A=Radius.*Radius*pi;

% Multiplier for volume calculation
CoMul=[Ledis(1)/2;(Ledis(1:end-1)+Ledis(2:end))/2;Ledis(end)/2];
Vol=CoMul'*A;

% Material prop
E=210e9*ones(numele,1);
Rho=7800*ones(numele,1);

% Calculation criterion
tol=1e-6;

%% 
[U, Lambda]=fem(E,Rho,Radius,Ledis,numele,eignum,fixeddofs,tol);
[BRadius,BLambda]=imodeshape(E,Rho,U,Ledis,numele,fixeddofs,Vol,tol);

%% Beam cross section plot
figure(1);clf;set(gca,'fontsize',32);
plot([0,L],[0,0],'-.k','Linewidth',2);hold on;
plot(x,Radius,'linewidth',4,'color','b');
h1=plot(x,-Radius,'linewidth',4,'color','b');
scatter(x,BRadius,10*length(BRadius),'or','MarkerFaceColor','r');
h2=scatter(x,-BRadius,10*length(BRadius),'or','MarkerFaceColor','r');
axis([0,L,-20,20]);
legend([h1,h2],'Given c/s Profile','Obtained c/s Profile');
xlabel('Axial Coordinate of the Beam (m)','fontsize',32);
ylabel('c/s of circular beam (m)','fontsize',32);

% 
%%
figure(2);axes('FontSize',32);
plot([0,L],[0,0],'-.k','Linewidth',2);hold on;
h3=plot(linspace(0,L,numele+1),U(1:2:end),'linewidth',3,'color',[1,0,0]);
grid on;
legend(h3,'Modeshape Provided');
axis([0,L,-2,2]);
xlabel('Axial Coordinate of the Beam (m)','fontsize',32);
ylabel('Normalized Modal (Axial) Displacement (m)','fontsize',32);