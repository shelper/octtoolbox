function LIMPAGO=GenerateLIMPAgo(CaliedData,  ZeroPadNum, Threshold)

AlineNum=size(CaliedData,2);
CaliedData=GetFreqComp(CaliedData',[AlineNum*Threshold,AlineNum/2+1], 2)';
LIMPAGO=abs(fft(CaliedData,ZeroPadNum));
LIMPAGO=LIMPAGO(2:end/2,:);


% LIMPAGOPos=GetFreqComp(ComplexImage',[AlineNum/10,AlineNum/2+1], 1)';
% LIMPAGONeg=GetFreqComp(ComplexImage',[AlineNum/2+1,AlineNum-AlineNum/10+2], 1)';
% 
% LIMPAGONeg=FiltImage(abs(LIMPAGONeg),LIMPFilterSize,LIMPFilterType)./IntensityImage;
% LIMPAGOPos=FiltImage(abs(LIMPAGOPos),LIMPFilterSize,LIMPFilterType)./IntensityImage;
% LIMPAGO=pi*(LIMPAGOPos-LIMPAGONeg);
