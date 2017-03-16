function [h]=plotbreadth(Breadth,L,scale,linewid)
	nele=length(Breadth);
	x=linspace(0,L,nele+1);
	
	xcoord=zeros(2*nele,1);
	xcoord(1:2:end)=x(1:nele);
	xcoord(2:2:end)=x(2:nele+1);
	
	breadth=zeros(2*nele,1);
	breadth(1:2:end)=Breadth/2;
	breadth(2:2:end)=Breadth/2;

	plot(xcoord,breadth,'linewidth',linewid,'color',[1,0,0]); hold on;
	h=plot(xcoord,-breadth,'linewidth',linewid,'color',[1,0,0]);
	maxbreadth=max(breadth);
	axis([0,L,-scale*maxbreadth,scale*maxbreadth]);
	grid on;
end
