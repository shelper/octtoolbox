function [MSSFlow, IntensityImage]=GenerateMSSFlow(CaliedData,FullRange, FilterSize, FilterType,ZeroPadNum)

% ComplexImage=ifft(CaliedData,ZeroPadNum);
% ComplexImage=FiltImage(ComplexImage, FilterSize, FilterType);
% IntensityImage = 0.25*abs(ComplexImage(:,3:end))+0.25*abs(ComplexImage(:,1:end-2))+0.5*abs(ComplexImage(:,2:end-1));
% % IntensityImage = sqrt(abs((ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)))));
% ComplexImage=diff(ComplexImage,1,2);
% MSSFlow = angle(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)));

ComplexImage=ifft(CaliedData,ZeroPadNum);
if ~FullRange
    ComplexImage=ComplexImage(2:end/2,:);
end

ComplexImage(:,2:end)=diff(ComplexImage,1,2);
MSSFlow  = angle(FiltImage(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)),FilterSize, FilterType));

IntensityImage = FiltImage(abs(ComplexImage(:,2:end))+abs(ComplexImage(:,1:end-1)),FilterSize, FilterType);
% IntensityImage=IntensityImage/median(max(IntensityImage))/10;

% ComplexImage=(ComplexImage, FilterSize, FilterType);
% PAMFlow = angle(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)));

