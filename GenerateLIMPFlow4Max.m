function LIMPFlow=GenerateLIMPFlow4Max(CaliedData, ZeroPadNum)

[PixelNum, AlineNum]=size(CaliedData); % get data size
ComplexImage=ifft(CaliedData,ZeroPadNum);% get the complex image by ifft

%% generate customized Linear Filter
LinearFilter=zeros(1,AlineNum);
LinearFilter(1:end/2+1)=0:2/AlineNum:1;
LinearFilter=repmat(LinearFilter, PixelNum,1);

%% filtering ComplexImage in lateral direction by the linear filter
LIMPFlowMap=ifft(fft(ComplexImage,[],2).*LinearFilter,[],2);

%% specify the filter to smooth the original amplitude and the final results
% change the filter size and type as you wish
ImageFilter=fspecial('average',3);
% use half of the complexImage to get the amplitude/intensityImage
IntensityImage=imfilter(abs(ComplexImage(2:end/2,:)),ImageFilter);

%% Flow reconstruction by ratio between the linearly filtered image (FlowMap) and the original amplitude image
% use the upper half as the pos part to map the positive flows and filtered
Pos=imfilter(abs(LIMPFlowMap(2:end/2,:)),ImageFilter)./IntensityImage;
% use the lower half as the neg part to map the negtive flows and filtered
Neg=imfilter(abs(LIMPFlowMap(end:-1:end/2+2,:)),ImageFilter)./IntensityImage;

%% get bidirectional flows  
LIMPFlow=pi*(Pos-Neg); % the flow rate is ranged from -pi to pi
figure;imshow(LIMPFlow,[-pi,pi]);colormap jet;