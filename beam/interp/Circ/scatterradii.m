function [h]=scatterradii(Radius,L,markersize)
	nele=length(Radius)-1;
	x=linspace(0,L,nele+1);
	
	plot(x,Radius,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	h=plot(x,-Radius,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize);
	grid on;
end
