function [A,Lambda]=PowerIter(psK,psM,K,M,Pow,Comul,Vol,Tol)
num=size(psK,1);

Mat=psK\psM;

A=ones(num,1);

Res=1.+Tol;
Change=1.0+Tol;

while(Res>Tol || Change>Tol)
    oldA=A;

    A=Mat*A;

    A=(abs(A).^(1/Pow));
    A=A/max(A);
    
    Lambda=(A'*(A.^Pow))/(A'*Mat*A);

    Res=norm(K*(A.^Pow)/Lambda-M*A)
    Change=norm(oldA-A);
    plot(A);hold on;pause;
       
end
obVol=Comul'*A;
rati=Vol/obVol;
A=rati*A;
Lambda=Lambda*(rati^(Pow-1));
end


