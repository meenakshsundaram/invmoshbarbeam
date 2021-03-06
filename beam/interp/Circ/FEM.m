function [U, Lambda]=FEM(E,Rho,Radius,L,numele,eignum,fixeddofs,tol)

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
KX=zeros(16*numele,1);KY=zeros(16*numele,1);KZ=zeros(16*numele,1);
MX=zeros(16*numele,1);MY=zeros(16*numele,1);MZ=zeros(16*numele,1);

ntriplets=0;
MI=pi*Radius.^4/4;
A=pi*Radius.^2;
for ele=1:numele
	eledofs=[2*ele-1 2*ele 2*ele+1 2*(ele+1)];
	gdofs=[ele ele+1];
	Keloc=Ke1*MI(gdofs(1))+Ke2*MI(gdofs(2));
	Keloc=Keloc*E(ele);
	Meloc=Me1*A(gdofs(1))+Me2*A(gdofs(2));
	Meloc=Meloc*Rho(ele);
    for i = 1:4
        for j=1:4
            ntriplets=ntriplets+1;
            KX(ntriplets)=eledofs(i); KY(ntriplets)=eledofs(j); KZ(ntriplets)=Keloc(i,j);
            MX(ntriplets)=eledofs(i); MY(ntriplets)=eledofs(j); MZ(ntriplets)=Meloc(i,j);
        end
    end
end

K=sparse(KX,KY,KZ,totdofs,totdofs);
M=sparse(MX,MY,MZ,totdofs,totdofs);

alldofs=1:totdofs;
freedofs=setdiff(alldofs,fixeddofs);

%solving
options.tol=tol;
options.issym=1;
options.isreal=1;
options.disp=0;
U=zeros(totdofs,1);

[eigvec eigval]=eigs(K(freedofs,freedofs),M(freedofs,freedofs),eignum,0,options);

U(freedofs)=eigvec(:,1);
Lambda=eigval(1,1);

%normalisation
[val,ind]=max(abs(U));
U=U/U(ind);

end
