function setDefaultsImagePrint(xLimits,yLimits,sz)
set(gca,'XTick',xLimits(1):4000:xLimits(2));
set(gca,'YTick',yLimits(1):2000:yLimits(2));
set(gca,'FontSize',sz);
axis equal

xlim(xLimits)
ylim(yLimits)
xticks(xLimits(1):4000:xLimits(2))
xticklabels({'418000','422000','426000','430000','434000'})
yticks(4514000:4000:4526000)
yticklabels({'4514000','4518000','4522000','4526000'})
xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',sz)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',sz)




