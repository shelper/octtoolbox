function GenerateHRatio(FileName,InitialHRatio, RawData,PrePar.ReferenceMode)
CaliedData = GetOCTData(FileName, 'ssOCTCam');
[PixelNum, AlineNum]=size(CaliedData);
CaliedData=RemoveRef(CaliedData,'mean');
% CaliedData=CaliedData.*repmat(barthannwin(PixelNum)+0.001,1,AlineNum);

%% phaase jittering correction
load(['C:\Users\zyuan\My Dropbox\research\OCTToolbox\FS\','0926.mat']);
[CaliedData,foo]=CorrectPixelShift(CaliedData,2,FSPhaseShift*0.64);

H0=GetFreqComp(CaliedData',[1, AlineNum/6],2)';
H1=CaliedData.*repmat([1,-1],[PixelNum,AlineNum/2]);%figure;plot(abs(fft(NewData,[],2))')
H1=GetFreqComp(H1',[1, AlineNum/6],2)';  
HRatio=mean(abs(H1),2)./mean(abs(H0),2);
figure;plot(HRatio);
fid = fopen(
FRHRatio=polyval(polyfit(1:1376,HRatio,5),1:1376);



