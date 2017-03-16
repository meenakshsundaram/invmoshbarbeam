clc;
clear;

numele=10;
eignum=1;
fixeddofs=[1,2];
L=1;

x=linspace(0,L,2*numele+1)';
x=x(2:2:end);

Breadth=sin(10*x)*5+10;
Depth=x+4;
Vol=sum(Breadth.*Depth*L/numele);

E=210e9*ones(size(x));
Rho=7800*ones(size(x));

flag=2;
tol=1e-9;
scale = 3;
linewid=2;

[U, Lambda]=FEM(E,Rho,Breadth,Depth,L,numele,eignum,fixeddofs,tol);
[BBreadth,BLambda]=IModeshape(E,Rho,U,Depth,L,numele,fixeddofs,Vol,tol*1e6,flag);

figure(1);axes('FontSize',16);
plot([0,L],[0,0],'-.k','Linewidth',linewid);hold on;

h1=plotbreadth(Breadth,L,scale,linewid);
h2=scatterbreadth(BBreadth,L,10);

legend([h1 h2],'Given c/s Profile','Obtained c/s Profile','location','best');

xlabel('Axial Coordinate of the Beam (m)','fontsize',16); 
ylabel('Side Length (m)','fontsize',16);

figure(2);axes('FontSize',16);
[h]=plotmodeshape(U,L,numele,linewid);
legend(h,'Modeshape Provided','location','best');
xlabel('Axial Coordinate of the Beam (m)','fontsize',16);
ylabel('Normalized Modal  Displacement (m)','fontsize',16);

