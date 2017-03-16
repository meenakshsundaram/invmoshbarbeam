function [h]=plotmodeshape(U,L,numele,linewid)
x=linspace(0,L,numele+1);
plot([0,L],[0,0],'-.k','Linewidth',2);hold on;
splineset=spapi(augknt([0 x L],4,2),[x x],[U(1:2:2*(numele+1))' U(2:2:2*(numele+1))']);
val=fnplt(splineset);
val(2,:)=val(2,:)/max(abs(val(2,:))); 
h=plot(val(1,:),val(2,:),'-b','Linewidth',linewid);hold on;
axis([0,L,-2,2]);
grid on;
end
