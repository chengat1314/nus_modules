function SetFont(FontName, TitleSize, TextSize, LabelSize, AxesSize);
% function SetFont(FontName, TitleSize, TextSize, LabelSize, AxesSize);
%   e.g.:  SetFont('Palatino', 18, 16, 18, 15);

h1=get(gcf,'Children');
level1=size(h1);
for i=(1:level1)
	if strcmp(get(h1(i),'type'),'axes')
		set(h1(i),'FontName',FontName,'FontSize',[AxesSize]);
		hx = get(h1(i),'XLabel');
		set(hx,'FontName',FontName,'FontSize',[LabelSize]);
		hy = get(h1(i),'YLabel');
		set(hy,'FontName',FontName,'FontSize',[LabelSize]);
		ht = get(h1(i),'Title');
		set(ht,'FontName',FontName,'FontSize',[TitleSize]);
	end
	h2=get(h1(i),'Children');
	level2=size(h2);
	for j=(1:level2)
		if strcmp(get(h2(j),'type'),'text')
			set(h2(j),'FontName',FontName,'FontSize',[TextSize]);
			set(h2(j),'HorizontalAlignment','center');
		end
	end
end
