clc;
clear;

numele=1005;
eignum=10;
fixeddofs=[1,2];
L=1;

x=linspace(0,L,2*numele+1)';
x=x(2:2:end);
Breadth=sin(10*x)*5+10;
Breadth=Breadth/max(Breadth);
Depth=x+4;

E=ones(size(x));
Rho=ones(size(x));

flag=1;
tol=1e-9;
[U, Lambda]=FEM(E,Rho,Breadth,Depth,L,numele,eignum,fixeddofs,tol);
[BBreadth,BLambda]=IModeshape(E,Rho,U,Depth,L,numele,fixeddofs,tol,flag);
