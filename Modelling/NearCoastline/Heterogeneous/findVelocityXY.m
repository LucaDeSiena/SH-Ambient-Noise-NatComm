function [tt,vdX,vdY]       =...
    findVelocityXY(Time,DisplacementX,DisplacementY,TSampling)

correctionTerm              =   10;
Fs                          =   1/TSampling/correctionTerm;
Nf                          =   50; 
Fpass                       =   2.2; 
Fstop                       =   3.6;
AmpTimeX                    =   DisplacementX(1:correctionTerm:end);
AmpTimeY                    =   DisplacementY(1:correctionTerm:end);
Time1                       =   Time(1:correctionTerm:end);
d = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',Fs);

% fvtool(d,'MagnitudeDisplay','zero-phase','Fs',Fs)
delay                       =   mean(grpdelay(d));
tt                          =   Time1(1:end-delay);
vdriftX                     =   filter(d,AmpTimeX)/TSampling;
vdriftY                     =   filter(d,AmpTimeY)/TSampling;
vdX                         =   vdriftX;
vdY                         =   vdriftY;
vdX(1:delay)                =   [];
vdY(1:delay)                =   [];
tt(1:delay)                 =   [];
vdX(1:delay)                =   [];
vdY(1:delay)                =   [];

% [pkp,lcp] = findpeaks(AmpTime);
% zcp = zeros(size(lcp));
% 
% [pkm,lcm] = findpeaks(-AmpTime);
% zcm = zeros(size(lcm));
% 
% subplot(2,1,1)
% plot(Time1,AmpTime,Time1([lcp lcm]),[pkp -pkm],'or')
% xlabel('Time (s)')
% ylabel('Displacement (cm)')
% grid
% 
% subplot(2,1,2)
% plot(tt,vd,Time1([lcp lcm]),[zcp zcm],'or')
% xlabel('Time (s)')
% ylabel('Speed (cm/s)')
% grid