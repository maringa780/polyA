function ax=labelheatmap(cgo,labels,hfig)

if nargin<3
ax=plot(cgo) ;
else
ax=plot(cgo,hfig) ;
end

if nargin>1
hold on

% Set grid lines
xline(ax,ax.XTick+0.5,'k-','Alpha',1)
yline(ax,ax.YTick+0.5,'k-','Alpha',1)


for i= ax.XTick
text(ax,repmat(i,numel(ax.YTick),1), ax.YTick, labels(:,i), ...
    'VerticalAlignment', 'middle','HorizontalAlignment','Center');
end
hold off
end