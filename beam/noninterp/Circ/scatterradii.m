function [h]=scatterRadius(Radius,L,markersize)
	nele=length(Radius);
	x=linspace(0,L,2*nele+1);
	
	xcoord=x(2:2:end);
	
	plot(xcoord,Radius,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	h=plot(xcoord,-Radius,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	grid on;
end
