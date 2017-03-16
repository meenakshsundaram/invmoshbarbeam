function [h]=scatterbreadth(Breadth,L,markersize)
	nele=length(Breadth);
	x=linspace(0,L,2*nele+1);
	
	xcoord=x(2:2:end);
	Breadth=Breadth/2;	
	plot(xcoord,Breadth,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	h=plot(xcoord,-Breadth,'o','MarkerFaceColor',[0,1,0],'MarkerSize',markersize); hold on;
	grid on;
end
