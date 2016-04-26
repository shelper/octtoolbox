function THQFlow=GenerateTHQFlow(RawData,ZeroPadNum)

% get the complex Image for phase and amplitude calculation
ComplexImage=fft(RawData,ZeroPadNum);

% get the phase difference between neibouring Alines, be awared that the
% Aline number is reduced by 1 due to the difference
ComplexImage=ComplexImage(2:end/2,2:end).*conj(ComplexImage(2:end/2,1:end-1));
Filter=fspecial('average', [3,3]);
ComplexImage=imfilter(ComplexImage,Filter);
% filter the ComplexImage, be careful the filtering is done in complex
% domain, it would be nicer if the size of the filtering is configurable


% get the flow image
THQFlow = angle(ComplexImage);

% get the intensity image which will be used to generate a binary mask by a
% thresholde to remove the phase noise from low signals. we can also use
% the IntensityImage from OCTImgConv, but be careful about the image size
% and filtering match between the flow and the IntensityImage
IntensityImage = sqrt(abs(ComplexImage));

% Masking by IntensityImage, here the IntensityThreshold is set to 10, it
% would be better if this can be configurable
IntensityThreshold=10;
THQFlow=THQFlow.*(IntensityImage>IntensityThreshold);




