clc;
clear;

numele=41;
eignum=2;
fixeddofs=[1,2];%,2*(numele)+1,2*(numele+1)];
L=1;

x=linspace(0,L,numele+1)';
Radius=sin(10*x)*5+20;
Vol=trapz(pi*Radius.*Radius)*L/numele;

E=210e9*ones(size(x));
Rho=7800*ones(size(x));
tol=1e-9;
scale = 3;
linewid=2;

[U, Lambda]=FEM(E,Rho,Radius,L,numele,eignum,fixeddofs,tol);
[BRadius,BLambda]=IModeshape(E,Rho,U,L,numele,fixeddofs,Vol,tol*1e3);

figure(1);axes('FontSize',16);
plot([0,L],[0,0],'-.k','Linewidth',linewid);hold on;

h1=plotradii(Radius,L,scale,linewid);
h2=scatterradii(BRadius,L,10);


legend([h1 h2],'Given c/s Profile','Obtained c/s Profile','location','best');

xlabel('Axial Coordinate of the Beam (m)','fontsize',16); 
ylabel('Side Length (m)','fontsize',16);

figure(2);axes('FontSize',16);
[h]=plotmodeshape(U,L,numele,linewid);
legend(h,'Modeshape Provided','location','best');
xlabel('Axial Coordinate of the Beam (m)','fontsize',16);
ylabel('Normalized Modal  Displacement (m)','fontsize',16);

