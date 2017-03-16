function [h]=plotdepth(Depth,L,scale,linewid)
	nele=length(Depth);
	x=linspace(0,L,nele+1);
	
	xcoord=zeros(2*nele,1);
	xcoord(1:2:end)=x(1:nele);
	xcoord(2:2:end)=x(2:nele+1);
	
	depth=zeros(2*nele,1);
	depth(1:2:end)=Depth/2;
	depth(2:2:end)=Depth/2;

	plot(xcoord,depth,'linewidth',linewid,'color',[1,0,0]); hold on;
	h=plot(xcoord,-depth,'linewidth',linewid,'color',[1,0,0]);
	maxdepth=max(depth);
	axis([0,L,-scale*maxdepth,scale*maxdepth]);
	grid on;
end
