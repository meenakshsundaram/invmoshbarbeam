function [A,Lambda]=PowerIter(psK,psM,K,M,Pow,Tol)
num=size(psK,1);
Mat=psK\psM;
A=ones(num,1);
Change=1.0+Tol;
Res=1.+Tol;
while(Res>Tol || Change>Tol)
    Aold=A;
    A=Mat*A;
    [val,ind]=max(abs(A));
    A=(abs(A/A(ind))).^(1/Pow);
    Lambda=A'*(A.^Pow)/(A'*Mat*A);
    Res=norm(K*(A.^Pow)*(1/Lambda)-M*A,2);
    Change=norm(Aold-A);    
end
end
