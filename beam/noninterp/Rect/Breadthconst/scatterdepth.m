function [h]=scatterDepth(Depth,L,markersize)
	nele=length(Depth);
	x=linspace(0,L,2*nele+1);
	
	xcoord=x(2:2:end);
	Depth=Depth/2;	
	plot(xcoord,Depth,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	h=plot(xcoord,-Depth,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	grid on;
end
