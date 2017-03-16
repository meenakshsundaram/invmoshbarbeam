function [Breadth,Lambda]=IModeshape(E,Rho,U,Depth,L,numele,fixeddofs,tol,flag)

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
KX=zeros(4*numele,1);KY=zeros(4*numele,1);KZ=zeros(4*numele,1);
MX=zeros(4*numele,1);MY=zeros(4*numele,1);MZ=zeros(4*numele,1);

ntriplets=0;
for ele=1:numele    
	eledofs=[2*ele-1 2*ele 2*ele+1 2*(ele+1)];
   Kvec=E(ele)*Ke*((Depth(ele)^3)/12)*U(eledofs);
   Mvec=Rho(ele)*Me*Depth(ele)*U(eledofs);
   for i = 1:4
        ntriplets=ntriplets+1;
        KX(ntriplets)=eledofs(i);KY(ntriplets)=ele;KZ(ntriplets)=Kvec(i);
        MX(ntriplets)=eledofs(i);MY(ntriplets)=ele;MZ(ntriplets)=Mvec(i);        
    end
end

K=sparse(KX,KY,KZ,totdofs,numele);
M=sparse(MX,MY,MZ,totdofs,numele);

%Setting up Boundary Conditions
alldofs=1:totdofs;
freedofs=setdiff(alldofs,fixeddofs);

%Pseudo Inverse
psK=K(freedofs,:)'*K(freedofs,:);
psM=K(freedofs,:)'*M(freedofs,:);
if(flag==1)
    [Breadth,Lambda]=PowerIter(psK,psM,K(freedofs,:),M(freedofs,:),1,tol);
else
    [Breadth,Lambda]=EigAna(psK,psM,K(freedofs,:),M(freedofs,:),tol);
end

end
