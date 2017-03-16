function [A,Lambda]=EigAna(psK,psM,K,M,tol)

%Eigenvalue Analysis
[eigvec eigval]=eig(full(psK),full(psM),'chol');
numele=size(psK,1);
%PostProcessing
%Finding Valid Areas
Num=0;
for ele=1:numele
    %the values should be reald
    if(isreal(eigvec(:,ele)) && isreal(eigval(ele,ele)))
        [val loc]=max(abs(eigvec(:,ele)));
        %Normalization
        eigvec(:,ele)=eigvec(:,ele)/eigvec(loc,ele);
        flag=0;
            for j=1:numele
                %after normalization none of the quantities should be
                %negative
                if(eigvec(j,ele)<tol)
                    flag=-1;
                    break;
                end
            end
            if(flag==0)
                err=norm(K*eigvec(:,ele)/eigval(ele,ele)-M*eigvec(:,ele),2);
                %check the equation solving
                if(err<tol)
                    Num=Num+1;
                    A(:,Num)=eigvec(:,ele);
                    Lambda(Num)=eigval(ele,ele);
                end
            end
            
    end
end
if(Num==0)
    A=[];
    Lambda=[];
end
end