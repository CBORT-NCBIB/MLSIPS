function h = hist2D(xdata,ydata,xbins,ybins)

h = zeros(numel(xbins),numel(ybins));

xInd = ceil((xdata(:)-xbins(1))/(xbins(end)-xbins(1))*numel(xbins));
xInd(xInd<=0) = 1;
xInd(xInd>numel(xbins)) = numel(xbins);

yInd = ceil((ydata(:)-ybins(1))/(ybins(end)-ybins(1))*numel(ybins));
yInd(yInd<=0) = 1;
yInd(yInd>numel(ybins)) = numel(ybins);

for ind = 1:numel(xdata)
    h(xInd(ind),yInd(ind)) = h(xInd(ind),yInd(ind)) + 1;
end