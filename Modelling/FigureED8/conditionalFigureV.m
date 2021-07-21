function true       =...
    conditionalFigureV(X,R,Az,xLimits,yLimits,P,interp,chosenCMap,...
    coordsE,C,F,tresholdR,velocityLimits,polarizationLimits,sz,...
    titleFigure)

conditionalX                    =   (X(:,1)<xLimits(1) |...
    X(:,1)>xLimits(2) |...
    X(:,2)<yLimits(1) | X(:,2)>yLimits(2));
X(conditionalX,:)=[];

conditionalP1                   =   (P(:,1)<xLimits(1) |...
    P(:,1)>xLimits(2) |...
    P(:,2)<yLimits(1) | P(:,2)>yLimits(2));

P(conditionalP1,:)              =   [];
cVect                           =   X(:,4);

PVect                           =   P(:,R);
PAz                             =   P(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

[xi,yi]                         =...
    meshgrid(xLimits(1):interp:xLimits(2), yLimits(1):interp:yLimits(2));
zi                              =   griddata(X(:,1),X(:,2),cVect,xi,yi);

true                            =   figure('Name',...
    'Group Velocity vs Residual Length and Azimuth, 2 s vs 0.2-1 Hz',...
    'NumberTitle','off','Position',[10 10 1300 560]);
ax1                             =   axes;
pcolor(xi,yi,zi);
setDefaultsImagePrint(xLimits,yLimits,velocityLimits,ax1,chosenCMap,sz)

hold on
plot(coordsE(:,2),coordsE(:,3),'d','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[1 1 1],'MarkerSize',8,'LineWidth',1)

mapshow(C,'FaceAlpha',0);
mapshow(F,'Color', 'black')

ax2                             =   axes;

scatter(ax2,P(:,1), P(:,2), 300, PVect,'Filled',...
    'MarkerEdgeColor','k','LineWidth',2);

hold on
conditionP                      =   P(:,R)>tresholdR;

quiver(P(conditionP,1), P(conditionP,2),...
    cX(conditionP),sX(conditionP),.1,'-w','LineWidth',3,...
    'ShowArrowHead','off')
quiver(P(conditionP,1), P(conditionP,2),...
    -cX(conditionP),-sX(conditionP),.1,'-w','LineWidth',3,...
    'ShowArrowHead','off')
setDefaultColorbarsPrint(ax1,ax2,polarizationLimits,chosenCMap,...
    titleFigure,sz)
xlim(xLimits)
ylim(yLimits)
hold off