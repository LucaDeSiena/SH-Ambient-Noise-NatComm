function RFigure(P,xLimits,yLimits,polarizationLimits,chosenCMap,sz,...
    coordsEarthquake)
plot(P(:,1),P(:,2),'o')
set(gca,'XTick',xLimits(1):4000:xLimits(2));
set(gca,'YTick',yLimits(1):2000:yLimits(2));
set(gca,'FontSize',12);

caxis(polarizationLimits)
xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',20)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',20)
colormap(chosenCMap)
axis equal
xlim(xLimits);
ylim(yLimits);
cb2                             =   colorbar;
cb2.Label.String                =   'R';
cb2.Label.FontSize              =   sz;
cb2.Label.FontWeight            =   'bold';
cb2.FontSize                    =   sz-2;
scatter(425800,4519200, 150, 'c','Filled',...
    'MarkerEdgeColor','w','LineWidth',1,'MarkerFaceColor','k');

if nargin == 7
    scatter(coordsEarthquake(1), coordsEarthquake(2), 250, 'c','Filled',...
    'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceColor','w');
end
    
