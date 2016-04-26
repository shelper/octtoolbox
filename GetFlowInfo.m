

function [T, SysT] = GetFlowInfo(PulseData)
Filter=ones(1,8)/8;
PulseDiff=imfilter(diff(PulseData),Filter);
[A,I]=max(PulseDiff);
tempPulseDiff=PulseDiff;
tempPulseDiff(I-5:I+5)=0;
[foo,tempI]=max(tempPulseDiff);
T=abs(tempI-I);

% get the right I before doing the following
B = 15;
tempData = PulseDiff(I-B:I+B);
C = 5;
D = mean(tempData([1:6,end-5:end]));
Params=fitgaussian(1:31, tempData,[A B C D]);

SysT= 3*Params(3);


