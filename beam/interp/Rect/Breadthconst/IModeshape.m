function [Depth,Lambda]=IModeshape(E,Rho,U,Breadth,L,numele,fixeddofs,tol)

Le=L/numele;
totdofs=2*(numele+1);
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

%Assembly
KX=zeros(8*numele,1);KY=zeros(8*numele,1);KZ=zeros(8*numele,1);
MX=zeros(8*numele,1);MY=zeros(8*numele,1);MZ=zeros(8*numele,1);

ntriplets=0;
for ele=1:numele    
	eledofs=[2*ele-1 2*ele 2*ele+1 2*(ele+1)];
	gdofs=[ele ele+1];
	Kvec1=E(ele)*Ke1*U(eledofs)*Breadth(gdofs(1))/12;
	Kvec2=E(ele)*Ke2*U(eledofs)*Breadth(gdofs(2))/12;
	Mvec1=Rho(ele)*Me1*U(eledofs)*Breadth(gdofs(1));
	Mvec2=Rho(ele)*Me2*U(eledofs)*Breadth(gdofs(2));
	for i = 1:4
        	ntriplets=ntriplets+1;
	        KX(ntriplets)=eledofs(i);KY(ntriplets)=gdofs(1);KZ(ntriplets)=Kvec1(i);
        	MX(ntriplets)=eledofs(i);MY(ntriplets)=gdofs(1);MZ(ntriplets)=Mvec1(i);
	end
	for i = 1:4
        	ntriplets=ntriplets+1;
	        KX(ntriplets)=eledofs(i);KY(ntriplets)=gdofs(2);KZ(ntriplets)=Kvec2(i);
        	MX(ntriplets)=eledofs(i);MY(ntriplets)=gdofs(2);MZ(ntriplets)=Mvec2(i);
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
[Depth,Lambda]=PowerIter(psK,psM,K(freedofs,:),M(freedofs,:),3,tol);
end
