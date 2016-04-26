function [TotalT, SysT, DiaT] = GetPulseInfo(PulseData, PulseSize)
PulseDiff=diff(PulseData);

Filter=ones(1,7)/7;
PulseDiff=imfilter(PulseDiff,Filter);


[foo,I1]=min(PulseDiff);
PulseDiff(max(1, I1-5):min(180, I1+5))=0;
[foo,I2]=min(PulseDiff);

if I1>I2
    foo = I1;
    I1 = I2;
    I2 = foo;
end
TotalT=I2-I1;

[foo,T1]=max(PulseData(I1:(I2+I1)/2));
[foo,T2]=min(PulseData((I2+I1)/2:I2));

DiaT=T2+round(TotalT/2)-T1;
SysT=TotalT-DiaT;



% % get the right I before doing the following
% B = 15;
% tempData = PulseDiff(I-B:I+B);
% C = 5;
% D = mean(tempData([1:6,end-5:end]));
% Params=fitgaussian(1:31, tempData,[A B C D]);
% 
% SysT= 4.47*Params(3);


