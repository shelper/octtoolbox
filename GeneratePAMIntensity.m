function IntensityImage=GeneratePAMIntensity(CaliedData,FilterSize, FilterType,ZeroPadNum)
ComplexImage=ifft(CaliedData,ZeroPadNum);
ComplexImage=FiltImage(ComplexImage, FilterSize, FilterType);
IntensityImage = abs((ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1))));
