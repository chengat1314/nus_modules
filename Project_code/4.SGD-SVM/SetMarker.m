function SetMarker(MarkerSize);
% function function SetMarker(MarkerSize);
%   e.g.:  SetMarker(4);

h1=get(gcf,'Children');
level1=length(h1);
for i=(1:level1)
	h2=get(h1(i),'Children');
	level2=length(h2);
	for j=(1:level2)
		if all(strcmp(get(h2(j),'type'),'line'))
			set(h2(j),'MarkerSize',[MarkerSize]);
		end
	end
end

