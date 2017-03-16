clc;
clear;

numele=1801;
eignum=10;
fixeddofs=[1,2];
L=1;

x=linspace(0,L,numele+1)';
Breadth=sin(10*x)*5+10;
Breadth=Breadth/max(Breadth);
Depth=x+4;
Depth=Depth/max(Depth);

E=210e9*ones(size(x));
Rho=7800*ones(size(x));
tol=1e-8;

[U, Lambda]=FEM(E,Rho,Breadth,Depth,L,numele,eignum,fixeddofs,tol);
[BDepth,BLambda]=IModeshape(E,Rho,U,Breadth,L,numele,fixeddofs,tol);
