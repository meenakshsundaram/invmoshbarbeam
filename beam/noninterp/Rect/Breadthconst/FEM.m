function [U, Lambda]=FEM(E,Rho,Breadth,Depth,L,numele,eignum,fixeddofs,tol)

Le=L/numele;
totdofs=2*(numele+1);

%Mass Matrix
Me=[	(13*Le)/35	(11*Le^2)/210	(9*Le)/70	-((13*Le^2)/420);
	(11*Le^2)/210	Le^3/105	(13*Le^2)/420	-(Le^3/140);
	(9*Le)/70	(13*Le^2)/420	(13*Le)/35	-((11*Le^2)/210);
	-((13*Le^2)/420)	-(Le^3/140)	-((11*Le^2)/210)	Le^3/105;];

%Stiffness Matrix
Ke=[	12/Le^3		6/Le^2		-(12/Le^3)	6/Le^2;
	6/Le^2		4/Le		-(6/Le^2)	2/Le;
	-(12/Le^3)	-(6/Le^2)	12/Le^3		-(6/Le^2);
	6/Le^2		2/Le		-(6/Le^2)	4/Le;	];

%Assembly
KX=zeros(16*numele,1);KY=zeros(16*numele,1);KZ=zeros(16*numele,1);
MX=zeros(16*numele,1);MY=zeros(16*numele,1);MZ=zeros(16*numele,1);

ntriplets=0;
for ele=1:numele
	eledofs=[2*ele-1 2*ele 2*ele+1 2*(ele+1)];
	Keloc=Ke*E(ele)*Breadth(ele)*(Depth(ele)^3)/12;
	Meloc=Me*Rho(ele)*Breadth(ele)*Depth(ele);
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
