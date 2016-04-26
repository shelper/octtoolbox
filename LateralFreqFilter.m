function FilteredImage=LateralFreqFilter(Image,FilterThresh,IntensityThresh)


fftImg=fft(Image,[],2);
i=int16(FilterThresh*size(Image,2))/2;
fftImg(:,end/2-i:end/2+i+2)=0;
fftImg=fftImg.*(abs(fftImg)>IntensityThresh);
FilteredImage=abs(ifft(fftImg,[],2));
% figure;imshow(20*log10(FiltedAbsImg),get_display_range(20*log10(FiltedAbsImg),[0.5,0.005]));
% figure;imshow((FiltedAbsImg),get_display_range((FiltedAbsImg),[0.5,0.005]));

% FftAbsImg=fft(AbsImg')';
% % ThreshedFftAbsImg=FftAbsImg.*(abs(FftAbsImg)>250);
% ThreshedFftAbsImg=FftAbsImg.*(abs(FftAbsImg)>3000);
% ThreshedAbsImg=abs(ifft(ThreshedFftAbsImg')');
% figure;imshow(20*log10(ThreshedAbsImg),get_display_range(20*log10(ThreshedAbsImg),[0.5,0.005]));
% % figure;imshow((ThreshedAbsImg),get_display_range((ThreshedAbsImg),[0.5,0.005]));
% 
% complex_image(1026:2048,:)=conj(complex_image(1024:-1:2,:));
% complex_image(1025,:)=0;
% complex_image(1,:)=0;
% RawData=ifft(complex_image);
% FftRawData=fft(RawData')';
% FftRawData(:,200:826)=0;
% % RawFiltedAbsImg=20*log10(abs(fft(ifft(FftRawData')',2048)));
% RawFiltedAbsImg=(abs(fft(ifft(FftRawData')')));
% RawFiltedAbsImg=RawFiltedAbsImg(2:1024,:);
% figure;imshow(20*log10(RawFiltedAbsImg),get_display_range(20*log10(RawFiltedAbsImg),[0.5,0.005]));
% % figure;imshow((RawFiltedAbsImg),get_display_range((RawFiltedAbsImg),[0.5,0.005]));
% 