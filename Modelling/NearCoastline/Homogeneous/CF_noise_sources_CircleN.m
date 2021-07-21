%% Anisotropic, viscoelastic SH-wave propagation at Campi Flegrei
% CODE for the simulation of SH waves in 2D for an anisotropic viscoelastic 
% medium. Reference publications: <https://library.seg.org/doi/abs/10.1190/1.1443792 
% Carcione and Cavallini 1996> + the staggered grid approach of <https://library.seg.org/doi/abs/10.1190/1.1444692 
% Carcione 1999>. Obviously, refer to Carcione (2014) also known as: Carcione, 
% J. M. (2007). Wave fields in real media: Wave propagation in anisotropic, anelastic, 
% porous and electromagnetic media: Handbook of geophysical exploration, seismic 
% exploration. _Handbook of Geophysical Exploration: Seismic Exploration_,  _38_.
% 
% Translated from the Fortran code.
%%
% 
%   Application: CAMPI FLEGREI
%   Area: 20x20 km
%   Dominant Frequency: 1 Hz
%   Filtered Frequency: 1 Hz
%   Frequency: 1 Hz
%   dx=V/dominant frequency: 2*10^3 (m/s)/ 0.7 (Hz) > 2800 -> 40 m
%   ddy=30000 m / 50 m = 600 nodes
%   ddx=30000 m / 50 m = 600 nodes
%   dt=(2/pi)* 50 (m)/ (3*10^3) (m/s) > 0.01 s -> 0.001
%

clear
close all
clc
warning off all
filename                        =   'coord_stations.xls';
delimiterIn                     =   ' ';
headerlinesIn                   =   0;
StationsNameCoords              =...
    importdata(filename,delimiterIn,headerlinesIn);
Coordslatlong                   =   StationsNameCoords.data;
NameStAll                       =   StationsNameCoords.textdata;
%%
% Location of the station on map
[Er,Nr,~]                       =   deg2utm(Coordslatlong(:,1),...
    Coordslatlong(:,2));
%% INPUTS
% Choose between isotropic (1), effective anisotropic (2) or scattering simulations 
% (3). In the last case choose the geometry between rectangular, circle or results.

chooseModel                         =   1;
chosenStation                       =   1;

% Name of the figures, will be added for each file;
labelFigures                        =   'CF_noisePolarization_CircleN_Homogeneous';

%% 
% Define velocity fluctuations:

rmsVelocityFluctuations             =   0;
correlationLength                   =   2;
%% 
% Define attenuation mechanisms:

Q2Average                           =   30;
Q4Average                           =   30;
%% Starting material properties
% Setting up velocities. You could follow <https://www.sciencedirect.com/science/article/abs/pii/S0377027319304044 
% Heap et al. 2020> if setting from rock physics e.g., V0=sqrt(5.4*10^9/2500)/1000=1.47 
% km/s. The ideal is to relate it to actual stiffness parameters so transform 
% them always - easy in the isotropic case. In this case, the frequency is too 
% low for those values .
% 
% The following will produce an anisotropic effective wavefield:
%% 
% * V0_WE = sqrt(9.5*10^9/rhoAverage)/1000; 
% * V0_SN = sqrt(12.5*10^9/rhoAverage)/1000; 
% * c460 = -1.5e+09;

rhoAverage                          =   2500;
rhoScattering                       =   1500;    
rhoSource                           =   2500;

V0_WE                               =   sqrt(2.25*10^10/rhoAverage)/1000;
V0_SN                               =   sqrt(2.25*10^10/rhoAverage)/1000;
V0                                  =   (V0_SN+V0_WE)/2;
c460                                =   0;
%% 
% In case you set up the scattering area:

if chooseModel >  2
    resultScattering                =   1;
    rectangleScattering             =   [];
    circleScattering                =   [];
    V0_SNScattering                 =   sqrt(3.25*10^10/rhoScattering)/1000;
    V0_WEScattering                 =   sqrt(3.25*10^10/rhoScattering)/1000;
    Q2Scattering                    =   30;
    Q4Scattering                    =   30;
end
%% Time/frequency inputs
% See calculations at the start. They are all in seconds - this simulation is 
% for *8* s. Start from the origin time in your SAC files (tData). For the frequency 
% inputs see calculations at the start - they are all in Hz. The relaxation frequency 
% and time define where you expect absorption to happen: here, set to *1* Hz so 
% that there is an effect.

dt                                  =   .001;
numberStep                          =   100000;
tData                               =   0;
Time                                =   tData+0:dt:tData+(numberStep-1)*dt;
fSource                             =   0.7;
fRelaxationTimes                    =   0.7;
relaxationTime                      =   1/(2*pi*fRelaxationTimes);
%% 
% Visualizes every nFrame steps

nFrame                              =   100;
%% Spatial inputs
% Set up the model nodes. Give the origin in UTM WGS84 as well as the number 
% of nodes and vectors in the 2D space. Stagger the grid with stresses (and all 
% other parameters) aginst displacements. 

dx                                  =   40;
dy                                  =   40;
origin                              =   [411100 4506100];
spatialSampling                     =   1/100;
nodesX                              =   750;
nodesY                              =   750;
%% Coordinates of the receivers, source and calculation of energy

% Location of the stations on grid and how big their markers:
rx                                  =   floor(Er-origin(1))/dx;
rz                                  =   floor(Nr-origin(2))/dy;
lR                                  =   200;

%% Build up space and variables with staggered vectorization

x                                   =   dx:dx:nodesX*dx;
X                                   =   origin(1) + x;
y                                   =   dy:dy:nodesY*dy;
Y                                   =   origin(2) + y;
xe                                  =   [x-dx/2,x(end)+dx/2];
ye                                  =   [y-dy/2,y(end)+dy/2];
Xe                                  =   origin(1) + xe;
Ye                                  =   origin(2) + ye;
[XX,YY]                             =   ndgrid(X,Y);
[XXx, YYx]                          =   ndgrid(Xe,Y);
[XXy, YYy]                          =   ndgrid(X,Ye);
%% 
% Field variable to update displacement (ux, uy, u2x, u2y):
%% 
% * e4 and e6 strain components;
% * s4 and s6 stress components;
% * e23 and e12 memory variables

u2x                                 =   zeros(nodesX+1,nodesY);
u2y                                 =   zeros(nodesX,nodesY+1);
ux                                  =   zeros(nodesX+1,nodesY);
uy                                  =   zeros(nodesX,nodesY+1);
s4                                  =   zeros(nodesX,nodesY);
s6                                  =   zeros(nodesX,nodesY);
e23                                 =   zeros(nodesX,nodesY);
e12                                 =   zeros(nodesX,nodesY);
rho                                 =   rhoAverage*ones(size(s4));
%% Scattering zones
% Set anisotropic zone for case 3. Create the volume you want to set as scattering 
% by hand. Also add the Rim, setting its thickness.

if chooseModel >  2
    
    if isempty(rectangleScattering) == 0
        scatterInitialX             =   nodesX/2 - 1000/dx;
        scatterFinishX              =   nodesX/2 + 1000/dx;
        scatterLX                   =   scatterInitialX:scatterFinishX;
        scatterLengthX              =   length(scatterLX);
        
        scatterInitialY             =   nodesY/2 - 5000/dy;
        scatterFinishY              =   nodesY/2 + 5000/dy;
        scatterLY                   =   scatterInitialY:scatterFinishY;
        scatterLengthY              =   length(scatterLY);
        
    end
    
    if isempty(circleScattering) == 0
        centerCaldera               =   [426500 4519500];
        rimThickness                =   [4000 6000];
        
        X1                          =   XX-centerCaldera(1);
        Y1                          =   YY-centerCaldera(2);
        R                           =   sqrt(X1.^2+Y1.^2);
        R1                          =...
            R>rimThickness(1) & R<rimThickness(2);
    end
    
    if isempty(resultScattering) == 0
        load anisotropicMesh3.mat
        mFint                       =   interp2(xq,yq,rFint,XX,YY);
        
        % In case the limits are outside of the grid interpolate better
        if find(isnan(mFint))
            mFint                   =   inpaintn(mFint);
        end
        maximumR                    =   0.3;
    end
end
%% Definition of the velocity fluctuations in space

csi                          =  fluctuationBuilder(s4,spatialSampling,...
    rmsVelocityFluctuations,correlationLength);
%% 
% Need to deal with both isotropic and anisotropic cases. As we use ndgrid the 
% program switches the input WE and SN!

if chooseModel ==  1
    csi_h                           =   V0*(1+csi);
elseif chooseModel ==  2
    csi_SN                          =   V0_WE*(1+csi);
    csi_WE                          =   V0_SN*(1+csi);
else
    csi_SN                          =   V0_WE*(1+csi);
    csi_WE                          =   V0_SN*(1+csi);
    if isempty(rectangleScattering) == 0
        csi_SN(scatterLY,scatterLX) =...
            V0_WEScattering*(1+csi(scatterLY,scatterLX));
        csi_WE(scatterLY,scatterLX)                  =...
            V0_SNScattering*(1+csi(scatterLY,scatterLX));
    end
    if isempty(circleScattering) == 0
        csi_SN(R1)                  =   V0_WEScattering*(1+csi(R1));
        csi_WE(R1)                  =   V0_SNScattering*(1+csi(R1));
    end
    if isempty(resultScattering) == 0
        xLimits=[418000 434000];
        yLimits=[4514000 4528000];
        x1=XX-xLimits(1);
        y1=YY-yLimits(1);
        z1=yLimits(2)-yLimits(1)+5000-0.9*x1;
        z2=yLimits(2)-yLimits(1)-2000-0.9*x1;
        z3=-3500+1*x1;
        z4=-2800+1*x1;

%             & XX>xLimits(1) & XX<xLimits(2) &...
%             YY>yLimits(1) & YY<yLimits(2));

%         upR                         =   (mFint>=maximumR &... 
%             XX>xLimits(1)+1000 & XX<xLimits(1)+ 1000 &...
%             YY>yLimits(2)-8000 & YY<yLimits(2)-2000);


        upR                         =   (mFint>=maximumR...
            & x1>0 & x1<24000 & y1>z2 & y1<z1) |...
            (x1>9000 & x1<10000 & y1>z3 & y1<z4);
%             (XX>=427000 & XX<=427700 & YY>=4520000 & YY<4521000);
        csi_SN(upR)                 =   V0_WEScattering*(1+csi(upR));
        csi_WE(upR)                 =   V0_SNScattering*(1+csi(upR));
    end
        
end
%% Definition of stiffness coefficients c44, c66,c46
% Standard anisotropic coefficients for the first two cases:

if chooseModel == 1
    
    c44                             =   rhoAverage.*(csi_h*1000).^2;
    c66                             =   rhoAverage.*(csi_h*1000).^2;
    c46                             =   zeros(size(s4));
     
elseif chooseModel==2
    
    c44                             =   rhoAverage.*(csi_WE*1000).^2;
    c66                             =   rhoAverage.*(csi_SN*1000).^2;
    c46                             =   c460*ones(nodesX,nodesY);

elseif chooseModel > 2
%% 
% Need to distinguish between background and high-scattering media; the standard 
% code had (for demonstration purposes):
%% 
% * c46Scattering = 0.5*sqrt(c44Scattering*c66Scattering);

    c44Average                      =   rhoAverage.*(csi_WE*1000).^2;
    c44Scattering                   =   rhoScattering.*(csi_WE*1000).^2;
    c66Average                      =   rhoAverage.*(csi_SN*1000).^2;
    c66Scattering                   =   rhoScattering.*(csi_SN*1000).^2;
    c46Average                      =   zeros(nodesX,nodesY);
    c46Scattering                   =   zeros(nodesX,nodesY);

end
%% Update field variables
% For the proper definition of each parameter look at tCarcione's book.

ts2Average                          =...
    (relaxationTime/Q2Average)*(sqrt(Q2Average*Q2Average+1)-1);
te2Average                          =...
    (relaxationTime/Q2Average)*(sqrt(Q2Average*Q2Average+1)+1);
phi2Average                         =   1/te2Average-1/ts2Average;
ts4Average                          =...
    (relaxationTime/Q4Average)*(sqrt(Q4Average*Q4Average+1)-1);
te4Average                          =...
    (relaxationTime/Q4Average)*(sqrt(Q4Average*Q4Average+1)+1);
phi4Average                         =   1/te4Average-1/ts4Average;

ts2                                 =   ts2Average*ones(nodesX,nodesY);
phi2                                =   phi2Average*ones(nodesX,nodesY);
ts4                                 =   ts4Average*ones(nodesX,nodesY);
phi4                                =   phi4Average*ones(nodesX,nodesY);

if chooseModel > 2
    ts2Scattering                       =...
        (relaxationTime/Q2Scattering)*(sqrt(Q2Scattering*Q2Scattering+1)-1);
    te2Scattering                       =...
        (relaxationTime/Q2Scattering)*(sqrt(Q2Scattering*Q2Scattering+1)+1);
    phi2Scattering                      =   1/te2Scattering-1/ts2Scattering;
    ts4Scattering                       =...
        (relaxationTime/Q4Scattering)*(sqrt(Q4Scattering*Q4Scattering+1)-1);
    te4Scattering                       =...
        (relaxationTime/Q4Scattering)*(sqrt(Q4Scattering*Q4Scattering+1)+1);
    phi4Scattering                      =   1/te4Scattering-1/ts4Scattering;
%% 
% Here we change the characteristics of the rectangle scattering area.

    if isempty(rectangleScattering) == 0
        
        c44                         =   c44Average;
        c44(scatterLX,scatterLY)=...
            c44Scattering(scatterLX,scatterLY);
        
        c66                         =   c66Average;
        c66(scatterLX,scatterLY)=...
            c66Scattering(scatterLX,scatterLY);
        
        c46                         =   c46Average;
        c46(scatterLX,scatterLY)=...
            c46Scattering(scatterLX,scatterLY);
        
        rho(scatterLX,scatterLY)=...
            rhoScattering*ones(scatterLengthX,scatterLengthY);
        
        ts2(scatterLX,scatterLY)=...
            ts2Scattering*ones(scatterLengthX,scatterLengthY);
        
        phi2(scatterLX,scatterLY)=...
            phi2Scattering*ones(scatterLengthX,scatterLengthY);
        
        ts4(scatterLX,scatterLY)=...
            ts4Scattering*ones(scatterLengthX,scatterLengthY);
        
        phi4(scatterLX,scatterLY)=...
            phi4Scattering*ones(scatterLengthX,scatterLengthY);
    end
%% 
% Here we change the characteristics of the rim.

    if isempty(circleScattering) == 0
        c44                         =   c44Average;
        c44(R1)                     =   c44Scattering(R1);
        c66                         =   c66Average;
        c66(R1)                     =   c66Scattering(R1);
        c46                         =   c46Average;
        c46(R1)                     =   c46Scattering(R1);
        ts2(R1)                     =   ts2Scattering;
        phi2(R1)                    =   phi2Scattering;
        ts4(R1)                     =   ts4Scattering;
        phi4(R1)                    =   phi4Scattering;
    
    end
%% 
% Here we change depending on the results. Using Petrosino & De Siena.

    if isempty(resultScattering) == 0
        c44                         =   c44Average;
        c44(upR)                    =   c44Scattering(upR);
        c66                         =   c66Average;
        c66(upR)                    =   c66Scattering(upR);
        c46                         =   c46Average;
        c46(upR)                    =   c46Scattering(upR);
        ts2(upR)                    =   ts2Scattering;
        phi2(upR)                   =   phi2Scattering;
        ts4(upR)                    =   ts4Scattering;
        phi4(upR)                   =   phi4Scattering;
    end

end

%% Setting up figure
% Choose the color map, size of markers and limits of visualization:

chosenCMap                          =   flipud(inferno);
sz                                  =   20;
xLimits=[418000 434000];
yLimits=[4514000 4528000];



% figure('Name',...
%     'Velocity fluctuations',...
%     'NumberTitle','off','Position',[10 10 2000 900]);
% hn                              =   pcolor(XX,YY,csi);
% setDefaultsImage(hn,xLimits,yLimits,chosenCMap,sz);
% 
% hold off
% nameFigure                          =...
%     cat(2,'CF_151007_Velocity_fluctuations_',...
%     num2str(rmsVelocityFluctuations),...
%     '_',num2str(correlationLength),'.tiff');
% print(nameFigure,'-dtiff','-r300');
%% 
% Here you need to insert the coordinates of the sources. You can 
% set it to a point (1), a circle around the area (2), or a vertical line
% (3).
optionSource                        =   2;

if optionSource == 1
    % Here you need to insert the coordinates of the source in degrees
    SN                              =   40 + 49.50/60;
    WE                              =   14 + 9.02/60;
    [Es,Ns,~]                       =   deg2utm(SN,WE);
    ix                              =   floor((Es-origin(1))/dx);
    iz                              =   floor((Ns-origin(2))/dy);

elseif optionSource == 2
    % Here you need to define the circumference giving centre and radius
    centerCaldera                   =   [426500 4521500];
    radiusSources                   =   [10000 10050];
    X1x                             =   XXx-centerCaldera(1);
    Y1x                             =   YYx-centerCaldera(2);
    X1y                             =   XXy-centerCaldera(1);
    Y1y                             =   YYy-centerCaldera(2);
    RSourcesx                       =   sqrt(X1x.^2+Y1x.^2);
    RSourcesy                       =   sqrt(X1y.^2+Y1y.^2);
    RSources1x                      =...
        RSourcesx>radiusSources(1) & RSourcesx<radiusSources(2);
    RSources1y                      =...
        RSourcesy>radiusSources(1) & RSourcesy<radiusSources(2);

elseif optionSource == 3
    % Here you need to insert the two x coords for the vertical line
        locationSource              =   [417000 417020];
        RSources1x                  =   XXx == locationSource(1);
        RSources1y                  =   XXy == locationSource(2);
end

%% 
% You could take the source intensity from many different approximations but 
% you always need to convert:
%% 
% * Velocity = Amp (counts) * AonD (volts)/counts / Gain * Sensitivity (volts/m*s).
%% 
% Energy is a function of Moment Magnitude:
%% 
% * log E = 5.24 + 1.44M -> E = exp(5.24 + 1.44*2.7)
%% 
% If 1.58997E-06 (V/counts) / (1 * 1500 (V/m/s)) then:
%% 
% * Amp = V/(AonD)*Gain*Sensitivity = V / 1.58997E-06* 1500.

M                                   =   0.01;
EGutenbergRichter                   =   10^(1.5*M + 4.8);
E                                   =   EGutenbergRichter;
differenceKmMeters                  =   1000;
Amp                                 =   sqrt(E);
rate                                =   numberStep/nFrame;
%% 
% The following is if you want it from the continuous wavelet transform:
%% 
% * [xrec,image,image2] =... sourceCWT('151007/151007-0910.CSOB.N.sac','Station_files/stazioni_151007.txt',... 
% fSource,50,54);
% * saveas(image,['CF_151007_',strcat(num2str(chooseModel)),'_Data_CSOB'])
% * saveas(image2,['CF_151007_',strcat(num2str(chooseModel)),'_CWT_CSOB'])
%% 
% As here I want to model SH waves from surface waves I take a morlet function:

dt2=dt/2;
[f,nw2]=morlet1(fSource,dt2,Amp,rhoSource);
nw2=floor(nw2);
nw=floor(nw2/2);

%% Show shear modulus variations

figure('Name',cat(2,labelFigures,'Shear Modulus'),...
    'NumberTitle','off','Position',[10 10 2000 900]);
hn                                  =   pcolor(XX,YY,c44);
setDefaultsImageShearModulus(hn,xLimits,yLimits,chosenCMap,sz);

for i=1:length(rx)
    rectangle('Position',[origin(1)+rx(i)*dx origin(2)+rz(i)*dy lR lR],...
        'FaceColor',[1 1 1])
    text(origin(1)+rx(i)*dx-2*lR,origin(2)+rz(i)*dy+1*lR,NameStAll{i},...
        'FontSize',14,'Color',[1 1 1],'FontWeight','bold','Rotation',-45)
end
hold off
nameFigure                          =...
    cat(2,cat(2,labelFigures,'Shear Modulus_'),...
    num2str(rmsVelocityFluctuations),...
    '_',num2str(correlationLength),'.tiff');
print(nameFigure,'-dtiff','-r300');

%% Set up wavefield and seismograms figures
% Set up amplitudes for seismograms, maximum number is 9:

AmpTimex                            =   zeros(numberStep,length(rx));
AmpTimey                            =   zeros(numberStep,length(rx));
%% 
% Set up the figure that shows the propagating wavefield:

figureWavefield                     =   figure('Name',...
    cat(2,labelFigures,'Wavefield evolution'),...
    'NumberTitle','off','Position',[1,200,900,800]);
panelWavefield                      =   surf(XXx,YYx,u2x/max(u2x(:)));
colormap(chosenCMap)
hold on
%% 
% Plotting for the scattering interface.

if chooseModel==3 && isempty(rectangleScattering) == 0
%% 
% Case of the rectangle:

    plot([origin(1)+scatterInitialX*dx origin(1)+scatterInitialX*dx],...
        [origin(2)+scatterInitialY*dy origin(2)+scatterFinishY*dy],'LineWidth',...
        2,'Color',[.5 .5 .5]);
    plot([origin(1)+scatterFinishX*dx origin(1)+scatterFinishX*dx],...
        [origin(2)+scatterInitialY*dy origin(2)+scatterFinishY*dy],'LineWidth',...
        2,'Color',[.5 .5 .5]);
    plot([origin(1)+scatterInitialX*dx origin(1)+scatterFinishX*dx],...
        [origin(2)+scatterInitialY*dy origin(2)+scatterInitialY*dy],'LineWidth',...
        2,'Color',[.5 .5 .5]);
    plot([origin(1)+scatterInitialX*dx origin(1)+scatterFinishX*dx],...
        [origin(2)+scatterFinishY*dy origin(2)+scatterFinishY*dy],'LineWidth',...
        2,'Color',[.5 .5 .5]);

elseif chooseModel==3 && isempty(circleScattering) == 0
%% 
% Case of the circular anomaly:

    viscircles([centerCaldera(1) centerCaldera(2)],...
        rimThickness(1),'Color',[1 1 1])
    viscircles([centerCaldera(1) centerCaldera(2)],...
        rimThickness(2),'Color',[1 1 1])

end
%% 
% Set up lighting, topography, source and stations:

light('Position',[0 0 0],'Style','infinite');
view(0,90)

C = shaperead('COSs.shp');
mapshow(C,'LineWidth',2,'FaceAlpha',0);

Faults= shaperead('FAULTS.shp');
mapshow(Faults,'Color', 'black','LineWidth',2)

for i=1:length(rx)
    rectangle('Position',[origin(1)+rx(i)*dx origin(2)+rz(i)*dy lR lR],...
        'FaceColor',[1 1 1])
    text(origin(1)+rx(i)*dx-2*lR,origin(2)+rz(i)*dy+1*lR,NameStAll{i},...
        'FontSize',14,'Color',[1 1 1],'FontWeight','bold','Rotation',-45)
end

xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',20)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',20)

hold off
%% 
% Set the data source for plotting. If you want to keep the colorbar constant 
% you can add:
%% 
% * caxis([-2*10^-4, 2*10^-4]);

set(panelWavefield,'xdatasource','XXx',...
    'ydatasource','YYx',...
    'zdatasource','u2x/max(u2x(:))')
set(panelWavefield,'FaceLighting','phong','FaceColor','interp',...
    'AmbientStrength',1,'SpecularStrength',1,...
    'DiffuseStrength',0.9,'SpecularColorReflectance',1,...
    'SpecularExponent',10,'FaceAlpha',1, 'linestyle','none')
setDefaultsImagePrint(xLimits,yLimits,sz)

c                                   =   colorbar;
c.Label.String                      =   'Normalized WE Displacement';
%% COMPUTATION
% Finite Differences weights:

weightX1                            =   9/(8*dx);
weightX2                            =   -1/(24*dx);
weightY1                            =   9/(8*dy);
weightY2                            =   -1/(24*dy);

% Absorbing boundaries parameters
r=1;
% r=0.999;
nab=140;
sab1(1:nab,1)=r.^(nab:-1:1);
% sab1(1:nab,1)=1;
% sab1(1:nab,1)=0;

% For u2x
sabtopX=repmat(sab1,1,nodesX);
sabbottomX=flipud(repmat(sab1,1,nodesX));
sableftX=repmat(sab1,1,nodesX+1)';
sabrightX=flipud(repmat(sab1,1,nodesX+1))';

% For u2y
sabtopY=repmat(sab1,1,nodesY+1);
sabbottomY=flipud(repmat(sab1,1,nodesY+1));
sableftY=repmat(sab1,1,nodesY)';
sabrightY=flipud(repmat(sab1,1,nodesY))';

%% 
% Starting time of real seismograms

t                                   =   0;
index                               =   0;
kx                                  =   3:nodesX-2;
kz                                  =   3:nodesY-2;
%% 
% Time stepping

tic
for n = 1:numberStep
    if mod(n*dt,10) == 0
        disp(num2str(n*dt))
    end
    t                               =   t+dt;
    % Absorbing boundaries
    
    % Horizontal stripes X
    u2x(1:nab,:)=u2x(1:nab,:).*sabtopX;
    u2x(nodesX-nab+1:nodesX,:)=u2x(nodesX-nab+1:nodesX,:).*sabbottomX;
    
    % Vertical stripes X
    u2x(:,1:nab)=u2x(:,1:nab).*sableftX;
    u2x(:,nodesY-nab+1:nodesY)=u2x(:,nodesY-nab+1:nodesY).*sabrightX;
    
    % Horizontal stripes Y
    u2y(1:nab,:)=u2y(1:nab,:).*sabtopY;
    u2y(nodesX-nab+1:nodesX,:)=u2y(nodesX-nab+1:nodesX,:).*sabbottomY;
    
    % Vertical stripes Y
    u2y(:,1:nab)=u2y(:,1:nab).*sableftY;
    u2y(:,nodesY-nab+1:nodesY)=u2y(:,nodesY-nab+1:nodesY).*sabrightY;
    
    %% 
% Strains

    e4                              =...
        weightY1*(u2y(kx,kz)-u2y(kx,kz-1)) +...
        weightY2*(u2y(kx,kz+1)-u2y(kx,kz-2));
    e6                              =...
        weightX1*(u2x(kx,kz)-u2x(kx-1,kz)) +...
        weightX2*(u2x(kx+1,kz)-u2x(kx-2,kz));
%% 
% Memory variables, absorption:

    f1                              =   2*ts2(kx,kz)-dt;
    f2                              =   2*ts2(kx,kz)+dt;
    ee                              =   e23(kx,kz);
    e23(kx,kz)                      =...
        (2*dt*ts2(kx,kz).*phi2(kx,kz).*e4+f1.*e23(kx,kz))./f2;
    e23(kx,kz)                      =   0.5*(e23(kx,kz)+ee);
    f1                              =   2*ts4(kx,kz)-dt;
    f2                              =   2*ts4(kx,kz)+dt;
    ee                              =   e12(kx,kz);
%% 
% Memory variables, anisotropy:

    e12(kx,kz)                      =...
        (2*dt*ts4(kx,kz).*phi4(kx,kz).*e6+f1.*e12(kx,kz))./f2;
    e12(kx,kz)                      =   0.5*(e12(kx,kz)+ee);
%% 
% Stresses:

    s4(kx,kz)                       =...
        c44(kx,kz).*(e4+e23(kx,kz))+c46(kx,kz).*e6;
    s6(kx,kz)                       =...
        c66(kx,kz).*(e6+e12(kx,kz))+c46(kx,kz).*e4;
%% 
% Differential stresses:

    ds4                             =...
        weightY1*(s4(kx,kz+1)-s4(kx,kz))+weightY2*(s4(kx,kz+2)-s4(kx,kz-1));
    ds6                             =...
        weightX1*(s6(kx+1,kz)-s6(kx,kz))+weightX2*(s6(kx+2,kz)-s6(kx-1,kz));
%% 
% Acceleration

    acceleration                    =   (ds4+ds6)./rho(kx,kz);
%% 
% Point source
if optionSource == 1
    if n<=nw
        n1                          =   n*2-1;
        source                      =   f(n1);
    elseif n>nw
        source                      =   0;
    end
% Extended source
elseif optionSource >= 2
    n1                              =   n*2-1;
    if n==1 || n1/nw2 == 1
        nSource                     =   1;
    else
        nSource                     =   mod(n1,nw2)+1;
    end
    source                          =   f(nSource);
end

%% 
% Euler equation:

    ux(kx,kz)                       =   2*u2x(kx,kz)-ux(kx,kz)+dt*dt*acceleration;
    ux(RSources1x)                   =   ux(RSources1x)+source;
    uy(kx,kz)                       =   2*u2y(kx,kz)-uy(kx,kz)+dt*dt*acceleration;
    uy(RSources1y)                   =   uy(RSources1y)+source;
%% 
% Update displacement

    uux                              =   u2x(kx,kz);
    uuz                              =   u2y(kx,kz);
    u2x(kx,kz)                       =   ux(kx,kz);
    u2y(kx,kz)                       =   uy(kx,kz);
    ux(kx,kz)                        =   uux;
    uy(kx,kz)                        =   uuz;
%% 
% Seismograms
    for i = 1:length(rx)
        AmpTimex(n,i)                    =...
            (u2x(floor(rx(i)+1),floor(rz(i)))...
            - u2x(floor(rx(i)-1),floor(rz(i))))/dx;
        AmpTimey(n,i)                    =...
            (u2y(floor(rx(i)),floor(rz(i)+1))...
            - u2y(floor(rx(i)),floor(rz(i)-1)))/dy;
    end
%% 
% Plot

    if isequal(mod(n,nFrame),0)
        index                       =   index+1;
        refreshdata(panelWavefield)
        drawnow
        F(index)                    =   getframe(figureWavefield);
        
    end
    
end
toc
%% 
% Set the figure where you plot the seismograms

seismogramsEvolution                =   figure('Name',...
    cat(2,labelFigures,'Model first seismogram'),'NumberTitle','off',...
    'Position',[1000,200,900,800]);

for j = 1:length(rx)
    [Time1,velocityWE(:,j),velocitySN(:,j)]       =   findVelocityXY(Time,AmpTimex(:,j),AmpTimey(:,j),dt);
end
title('First Seismogram')
subplot(2,1,1) 
plot(Time1,velocityWE(:,1),'Color',[0.64,0.08,0.18],'LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Normalized Velocity WE ','FontSize',12,'FontWeight','normal','Color','k')

subplot(2,1,2) 
plot(Time1,velocitySN(:,1),'Color',[0.64,0.08,0.18],'LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Normalized Velocity SN ','FontSize',12,'FontWeight','normal','Color','k')

labelVideo                      =   cat(2,labelFigures,'.avi');
labelMatFig                     =   cat(2,labelFigures,'_Modelled.fig');

%% 
% Create Video and Figures

video                               =...
    VideoWriter(labelVideo,'Motion JPEG AVI' );
video.FrameRate                     =   3;
open(video);
writeVideo(video,F);
close(video);
saveas(seismogramsEvolution,labelMatFig)

save(cat(2,labelFigures,'.mat'),'velocityWE','velocitySN');