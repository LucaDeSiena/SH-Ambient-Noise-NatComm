function csi_h                  =   fluctuationBuilder(u2,...
    spatialSampling,rmsVelocityFluctuations,correlationLength)

su2                             =   size(u2);
spatialWavenumberX              =   zeros(su2);
spatialWavenumberX(1,:)         =...
    spatialSampling:spatialSampling:(su2(2))/100;
spatialWavenumberX(:,1)         =...
    spatialSampling:spatialSampling:(su2(1))/100;

m1                              =   zeros(su2(1),1);
for i = 2:su2(2)
    m1(i)                       =   spatialSampling*i;
    spatialWavenumberX(2:end,i) =...
        sqrt(m1(i).^2+spatialWavenumberX(2:end,1).^2);
end
%% 
% The power spectral density function is exponential.

Pexp                            =...
    4*pi*rmsVelocityFluctuations^2*correlationLength^2./...
    (1+correlationLength^2*spatialWavenumberX.^2).^2;
phase                           =...
    2*pi.*rand(size(spatialWavenumberX));
ePhase                          =   (exp(1i*phase));
%% 
% Velocity fluctuations respect to the average.

csi_h                           =...
    1/4/pi^2*abs(fft2(sqrt(Pexp).*ePhase));
csi_h                           =   csi_h-mean(csi_h(:));