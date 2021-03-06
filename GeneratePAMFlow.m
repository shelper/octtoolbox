function [PAMFlow, IntensityImage]=GeneratePAMFlow(CaliedData,FullRange, FilterSize, FilterType,ZeroPadNum, Kasai)

ComplexImage=fft(CaliedData,ZeroPadNum);

%% filtering ComplexImage in lateral direction by the linear filter
% AlineNum=size(ComplexImage,2); % get data size
% ComplexImage=GetFreqComp(ComplexImage',[AlineNum/20, AlineNum/2],3)';
% Filter=zeros(1,AlineNum);
% Filter(1:end/2+1)=0:2/AlineNum:1;
% Filter=repmat(Filter, PixelNum,1);
% LIMPFlowMap=ifft(fft(ComplexImage,[],2).*Filter,[],2);


if ~FullRange
    ComplexImage=ComplexImage(2:end/2,:);
else
    ComplexImage=fftshift(ComplexImage,1);
end

% ComplexImage=sqrt(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)));
% ComplexImage=sqrt(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)));
% ComplexImage(:,2:end)=sqrt(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1))); 
% Angle=angle(ComplexImage);
% Angle= (Angle>0.5*pi)*(-pi)+(Angle<-0.5*pi)*pi;
% ComplexImage=ComplexImage.*exp(1i*Angle);

% Phase=GetFreqComp(angle(ComplexImage'),[AlineNum/2+1,AlineNum/2+1],2)'*6;
% ComplexImage=abs(ComplexImage).*exp(1i*Phase);
% 
% ComplexImage(:,2:end)=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));
% ComplexImage(:,2:end)=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));     
% ComplexImage(:,2:end)=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));     


if Kasai
    ComplexImage(:,2:end)=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));
    ComplexImage = FiltImage(ComplexImage,FilterSize, FilterType);
    PAMFlow = angle(ComplexImage);
else
%     ComplexImage(:,2:end)=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));
    PAMFlow = FiltImage(angle(ComplexImage),FilterSize, FilterType);
end

% figure;plot(PAMFlow(422,:))
IntensityImage = sqrt(abs(FiltImage(ComplexImage,FilterSize, FilterType)));
% IntensityImage=IntensityImage/median(max(IntensityImage))/10;

% ComplexImage=(ComplexImage, FilterSize, FilterType);
% PAMFlow = angle(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)));
% IntensityImage = sqrt(abs((ComplexImage(:,2:end)).*conj(ComplexImage(:,1:end-1))))/2;