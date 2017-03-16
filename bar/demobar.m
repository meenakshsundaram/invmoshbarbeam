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
clear;close all;

% Number of elements
numele=12;

% Element properties
E=210e9*ones(numele,1);
Rho=7800*ones(numele,1);

% Total length of bar
L=1;

% node position
x=linspace(0,L,numele+1)';

% area variation
A=5*sin(x*10)+10;

% integrate to find total volume
Vol=trapz(A)*L/numele;

% fixed degrees of freedom
fixeddofs=[1,2*numele+1];

% eigenmode desired
eignum=2;

% tolerance for problem
tol=1e-9;

%% Generate mode shape
[U,Lambda]=fem(E,Rho,A,L,numele,eignum,fixeddofs,tol);

%% Solve for the inverse mode shape
[BA,BLambda]=imodeshape(E,Rho,U,L,numele,fixeddofs,Vol,tol,1);

%% Plot
figure(1);clf;set(gca,'fontsize',32);
plot([0,L],[0,0],'-.k','Linewidth',2);hold on;
plot(x,sqrt(A),'linewidth',4,'color','b');
h1=plot(x,-sqrt(A),'linewidth',4,'color','b');
scatter(x,sqrt(BA),10*length(BA),'or','MarkerFaceColor','r');
h2=scatter(x,-sqrt(BA),10*length(BA),'or','MarkerFaceColor','r');
axis([0,1,-10,10]);
legend([h1,h2],'Given c/s Profile','Obtained c/s Profile');
xlabel('Axial Coordinate of the Bar (m)','fontsize',32);
ylabel('c/s of bar (m)','fontsize',32);

%%
figure(2);axes('FontSize',32);
plot([0,L],[0,0],'-.k','Linewidth',2);hold on;
h3=plot(linspace(0,L,2*numele+1),U,'linewidth',3,'color',[1,0,0]);
grid on;
legend(h3,'Modeshape Provided');
axis([0,L,-2,2]);
xlabel('Axial Coordinate of the Bar (m)','fontsize',32);
ylabel('Normalized Modal (Axial) Displacement (m)','fontsize',32);