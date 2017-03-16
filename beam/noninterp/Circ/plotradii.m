function [h]=plotradii(Radius,L,scale,linewid)
	nele=length(Radius);
	x=linspace(0,L,nele+1);
	
	xcoord=zeros(2*nele,1);
	xcoord(1:2:end)=x(1:nele);
	xcoord(2:2:end)=x(2:nele+1);
	
	radii=zeros(2*nele,1);
	radii(1:2:end)=Radius;
	radii(2:2:end)=Radius;

	plot(xcoord,radii,'linewidth',linewid,'color',[1,0,0]); hold on;
	h=plot(xcoord,-radii,'linewidth',linewid,'color',[1,0,0]);
	maxradii=max(radii);
	axis([0,L,-scale*maxradii,scale*maxradii]);
	grid on;
end
