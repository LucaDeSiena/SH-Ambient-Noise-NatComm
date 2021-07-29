%% Inputs
close all; clc; clear;
filename                        =   'coord_stations.xls';
delimiterIn                     =   ' ';
headerlinesIn                   =   0;
StationsNameCoords              =...
    importdata(filename,delimiterIn,headerlinesIn);
Coordslatlong                   =   StationsNameCoords.data;
NameStAll                       =   StationsNameCoords.textdata;

C                               =   shaperead('COSs.shp');
F                               =   shaperead('FAULTS.shp');
chosenCMap                      =   flipud(inferno);
% Inputs
xLimits                         =   [418000 434000];
yLimits                         =   [4514000 4528000];

% polarizationLimits              =   [0.2 0.6];
R                               =   3; %R polarization parameter
Az                              =   1; %Az polarization parameter
interp                          =   100;
sz                              =   24;
lR                              =   200;
dx                              =   40;
dy                              =   40;

tresholdR                       =   0.25;

%%
% Location of the station on map
[Er,Nr,~]                       =   deg2utm(Coordslatlong(:,1),...
    Coordslatlong(:,2));

%%
% Select stations where you want azimuths
[Obs,NameSt]                    =...
    importHrNoise2('2017_2009_stat.xls',NameStAll,Er,Nr);
AzObs                           =   Obs(:,Az)-90;
conditionP                      =   Obs(:,R)>tresholdR;
PObs                            =   Obs(conditionP,:);
ObsNameSt                       =   NameSt(conditionP,:);
%% Import coordinates stations
% Location of the station on map
origin                          =   [411100 4506100];
rx                              =   PObs(:,7);
rz                              =   PObs(:,8);
PObsAz                          =   PObs(:,Az)-90;

%% Interpolated map 0.2-1 - Line Homogeneous
% Reference earthquake
% Meshing
x                               =   xLimits(1):100:xLimits(2);
y                               =   yLimits(1):100:yLimits(2);

[XX,YY]                         =   ndgrid(x,y);
c44                             =   10^9*ones(size(XX));
AzR1                            =   load('AZ_Rlength_Line_Hom.dat');
labelFigures                    =   'CF_noisePolarization_Line_Homogeneous';

for i = 1:length(ObsNameSt)
    for j = 1:length(NameStAll)
        if isequal(ObsNameSt(i),NameStAll(j))
            AzR(i,:)            =   AzR1(j,:);
            break
        end
    end
end

PAz                             =   AzR(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

figure('Name',cat(2,labelFigures,'Shear Modulus with Azimuths'),...
    'NumberTitle','off','Position',[10 10 2000 900]);
hn                                  =   pcolor(XX,YY,c44);
setDefaultsImageShearModulus(hn,xLimits,yLimits,chosenCMap,sz);

for i=1:length(rx)
    rectangle('Position',[rx(i) rz(i) lR lR],...
        'FaceColor',[1 1 1])
    text(rx(i)-2*lR,rz(i)+1*lR,NameStAll{i},...
        'FontSize',14,'Color',[1 1 1],'FontWeight','bold','Rotation',-45)
end
hold on
quiver(rx+100, rz+100,cX,sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
quiver(rx+100, rz+100,-cX,-sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
hold off
nameFigure                          =...
    cat(2,cat(2,labelFigures,'Azimuths'),'.tiff');
print(nameFigure,'-dtiff','-r300');
AzimuthDifference                   =   PAz - PObsAz;
%% Interpolated map 0.2-1 - Line Continuous
load c44.mat
[nodesX, nodesY]=size(c44);
x                                   =   dx:dx:nodesX*dx;
X                                   =   origin(1) + x;
y                                   =   dy:dy:nodesY*dy;
Y                                   =   origin(2) + y;
[XX,YY]                             =   ndgrid(X,Y);

AzR1                            =   load('AZ_Rlength_Line_Con.dat');
labelFigures                    =   'CF_noisePolarization_Line_Continuous';

for i = 1:length(ObsNameSt)
    for j = 1:length(NameStAll)
        if isequal(ObsNameSt(i),NameStAll(j))
            AzR(i,:)            =   AzR1(j,:);
            break
        end
    end
end


PAz                             =   AzR(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

figure('Name',cat(2,labelFigures,'Shear Modulus with Azimuths'),...
    'NumberTitle','off','Position',[10 10 2000 900]);
hn                                  =   pcolor(XX,YY,c44);
setDefaultsImageShearModulus(hn,xLimits,yLimits,chosenCMap,sz);

for i=1:length(rx)
    rectangle('Position',[rx(i) rz(i) lR lR],...
        'FaceColor',[1 1 1])
    text(rx(i)-2*lR,rz(i)+1*lR,NameStAll{i},...
        'FontSize',14,'Color',[0 0 0],'FontWeight','bold','Rotation',-45)
end
hold on
quiver(rx+100, rz+100,cX,sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
quiver(rx+100, rz+100,-cX,-sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
hold off
nameFigure                          =...
    cat(2,cat(2,labelFigures,'Azimuths'),'.tiff');
print(nameFigure,'-dtiff','-r300');
AzimuthDifference(:,2)              =   PAz - PObsAz;

%% Interpolated map 0.2-1 - Circle North Homogeneous
% Reference earthquake
% Meshing
x                               =   xLimits(1):100:xLimits(2);
y                               =   yLimits(1):100:yLimits(2);

[XX,YY]                         =   ndgrid(x,y);
c44                             =   10^9*ones(size(XX));
AzR1                            =   load('AZ_Rlength_CircleN_Hom.dat');
labelFigures                    =   'CF_noisePolarization_CircleN_Homogeneous';

for i = 1:length(ObsNameSt)
    for j = 1:length(NameStAll)
        if isequal(ObsNameSt(i),NameStAll(j))
            AzR(i,:)            =   AzR1(j,:);
            break
        end
    end
end

PAz                             =   AzR(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

figure('Name',cat(2,labelFigures,'Shear Modulus with Azimuths'),...
    'NumberTitle','off','Position',[10 10 2000 900]);
hn                                  =   pcolor(XX,YY,c44);
setDefaultsImageShearModulus(hn,xLimits,yLimits,chosenCMap,sz);

for i=1:length(rx)
    rectangle('Position',[rx(i) rz(i) lR lR],...
        'FaceColor',[1 1 1])
    text(rx(i)-2*lR,rz(i)+1*lR,NameStAll{i},...
        'FontSize',14,'Color',[1 1 1],'FontWeight','bold','Rotation',-45)
end
hold on
quiver(rx+100, rz+100,cX,sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
quiver(rx+100, rz+100,-cX,-sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
hold off
nameFigure                          =...
    cat(2,cat(2,labelFigures,'Azimuths'),'.tiff');
print(nameFigure,'-dtiff','-r300');
AzimuthDifference(:,3)              =   PAz - PObsAz;

%% Interpolated map 0.2-1 - Circle North Continuous
load c44.mat
[nodesX, nodesY]=size(c44);
x                                   =   dx:dx:nodesX*dx;
X                                   =   origin(1) + x;
y                                   =   dy:dy:nodesY*dy;
Y                                   =   origin(2) + y;
[XX,YY]                             =   ndgrid(X,Y);

AzR1                            =   load('AZ_Rlength_CircleN_Con.dat');
labelFigures                    =   'CF_noisePolarization_CircleN_Continuous';

for i = 1:length(ObsNameSt)
    for j = 1:length(NameStAll)
        if isequal(ObsNameSt(i),NameStAll(j))
            AzR(i,:)            =   AzR1(j,:);
            break
        end
    end
end


PAz                             =   AzR(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

figure('Name',cat(2,labelFigures,'Shear Modulus with Azimuths'),...
    'NumberTitle','off','Position',[10 10 2000 900]);
hn                                  =   pcolor(XX,YY,c44);
setDefaultsImageShearModulus(hn,xLimits,yLimits,chosenCMap,sz);

for i=1:length(rx)
    rectangle('Position',[rx(i) rz(i) lR lR],...
        'FaceColor',[1 1 1])
    text(rx(i)-2*lR,rz(i)+1*lR,NameStAll{i},...
        'FontSize',14,'Color',[0 0 0],'FontWeight','bold','Rotation',-45)
end
hold on
quiver(rx+100, rz+100,cX,sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
quiver(rx+100, rz+100,-cX,-sX,.1,'-r','LineWidth',4,...
    'ShowArrowHead','off')
hold off
nameFigure                          =...
    cat(2,cat(2,labelFigures,'Azimuths'),'.tiff');
print(nameFigure,'-dtiff','-r300');
AzimuthDifference(:,4)              =   PAz - PObsAz;

sAz=sqrt(sum(AzimuthDifference.^2));