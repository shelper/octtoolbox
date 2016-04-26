function [AngleRange,VoltageRange]=CalculateVoltageRange(AlineNum,PassTime,Angle2Voltage,GalvoOffset)
%% calculate the mirror offset
% VoltageRange=CalculateVoltageRange(AlineNum,PassTime,Angle2Voltage,GalvoOffset)
% Angle2Voltage is 0.3 for Angle in degree unit
if nargin==2
    GalvoOffset=1500;
    Angle2Voltage=0.3;%Angle in degree;
end
Wavelength=1.06;
PasslengthStep=Wavelength/4/PassTime;

AngleRange=atan((AlineNum-1)*PasslengthStep/GalvoOffset/2);
AngleRange=AngleRange*180/pi;
VoltageRange=AngleRange*Angle2Voltage;


% AngleRange=10;
% Theta1=-AngleRange/2:AngleRange/(AlineNum-1):AngleRange/2;
% PathLength=2*GalvoOffset*(tan(pi/4+Theta1*pi/180)-tan(pi/4));
% 
% figure;plot(Theta1,PathLength);xlabel('GalvoAngle'); ylabel('PassLength');
% figure;plot(Theta1(2:end),diff(PathLength)/1.060*2);xlabel('GalvoAngle'); ylabel('PhaseChangePerAline (x pi)');
% 

%% this part calculates path change by the glass transmission 
% clear all;
% close all;
% %for transmission phase shift
% RefractiveIndex=1.7;
% GlassThickness=1000;
% i=1;
% for Theta1=-8:0.1:8
%     Theta1=Theta1*0.01745;
%     Theta2=asin(sin(Theta1)/RefractiveIndex);
%     PathLength(i)=RefractiveIndex*GlassThickness/cos(Theta2)-GlassThickness/cos(Theta2)*cos(Theta1-Theta2);
%     PosShift(i)=GlassThickness/cos(Theta2)*sin(Theta1-Theta2);
%     i=i+1;
% end
% figure;plot(PathLength);
% hold on; plot(PosShift,'r');
% 
% %for transmission phase shift
% clear all;
% close all;