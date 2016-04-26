function [LIMPFlow,IntensityForThresh]=GenerateAbsLIMPFlow(CaliedData, LIMPFilterSize,LIMPFilterType, ZeroPadnum, MSS)

IntensityForThresh=abs(FiltImage(ifft(CaliedData,ZeroPadnum),LIMPFilterSize,LIMPFilterType));
if MSS
    CaliedData=diff(CaliedData,1,2);
    IntensityForThresh=((IntensityForThresh(:,2:end))+(IntensityForThresh(:,1:end-1)))/2;
end
% IntensityForThresh=IntensityForThresh/median(max(IntensityForThresh))/10;


Intensity=(ifft(CaliedData,ZeroPadnum));
AlineNum=size(CaliedData,2);
fftData=fft(CaliedData,[],2);
Filter=zeros(AlineNum,1);Filter(2:end)=triang(AlineNum-1);
LIMPFlow=(ifft(ifft(fftData*sparse(diag(Filter)),[],2),ZeroPadnum));
% LIMPFlow=abs(FiltImage(LIMPFlow,  LIMPFilterSize,LIMPFilterType))./abs(FiltImage(Intensity, LIMPFilterSize,LIMPFilterType))*pi; 
LIMPFlow=abs(FiltImage(LIMPFlow./Intensity,  LIMPFilterSize,LIMPFilterType))*pi; 
LIMPFlow=min(LIMPFlow,pi);

