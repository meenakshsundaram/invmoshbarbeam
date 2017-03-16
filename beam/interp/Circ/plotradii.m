function [h]=plotradii(Radius,L,scale,linewid)
	nele=length(Radius)-1;
	x=linspace(0,L,nele+1);
	
	plot(x,Radius,'linewidth',linewid,'color',[1,0,0]); hold on;
	h=plot(x,-Radius,'linewidth',linewid,'color',[1,0,0]);
	maxRadii=max(Radius);
	axis([0,L,-scale*maxRadii,scale*maxRadii]);
	grid on;
end
