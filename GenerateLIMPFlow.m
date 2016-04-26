function [LIMPFlow,IntensityImage]=GenerateLIMPFlow(CaliedData, FullRange, LIMPFilterSize,LIMPFilterType, ZeroPadNum, MSS)

if MSS
    CaliedData(:,2:end)=diff(CaliedData,1,2);
end

% if ~FullRange
%     ComplexImage=ComplexImage(2:end/2,:);
% end
%% generate customized Linear Filter
% [ZeroPadNum, AlineNum]=size(CaliedData);
AlineNum = size(CaliedData,2);
% Filter=fliplr(fftshift(2/AlineNum-1:2/AlineNum:1));
Filter=zeros(1,AlineNum);
Filter(1:end/2+1)=0:2/AlineNum:1;
% Filter=Filter*1.2;
Filter=repmat(Filter, ZeroPadNum,1);

%% method 1, filtering interferogram data
% foo=ones(size(CaliedData));
% foo= mean(sum(abs(CaliedData)))*foo./repmat(sum(abs(CaliedData)), 2048,1);
% CaliedData=CaliedData.*foo;
ComplexImage=ifft(CaliedData,ZeroPadNum);
% figure;plot(sum(abs(CaliedData2)))
% CaliedData=CaliedData1;
% fftData=fft(CaliedData,[],2);
% LIMPFlowMap=ifft(ifft(fftData.*Filter,[],2),ZeroPadNum);

%% method 2, filtering ComplexImage
ComplexImage=exp(1i.*angle(ComplexImage));
fftData=fft(ComplexImage,[],2);
LIMPFlowMap=ifft(fftData.*Filter,[],2);

%% different filtering approaches.
IntensityImage=abs(ComplexImage(2:end/2,:));
% Neg=abs(LIMPFlowMap(2:end/2,:))./IntensityImage;
% Pos=abs(LIMPFlowMap(end:-1:end/2+2,:))./IntensityImage;

Filter=fspecial('average',1);
Neg=imfilter(abs(LIMPFlowMap(2:end/2,:)),Filter);
Pos=imfilter(abs(LIMPFlowMap(end:-1:end/2+2,:)),Filter);

% PosMask= (Pos>=(Neg./min((Neg+1)/2,1)));%in order to resolve the part close to PI
% NegMask= (Neg>(Pos./min((Pos+1)/2,1)));%in order to resolve the part close to PI
PosMask= (Pos>=Neg);%in order to resolve the part close to PI
NegMask= (Neg>Pos);%in order to resolve the part close to PI


% PosMask=Pos>Neg;
% NegMask=Neg>Pos;
% when LIMPPos and LIMPNeg close to 1 (high freq) then, PosMask


%% specific fast flow processing
% PosNegDiff=(Pos-Neg)./(Pos+Neg);
% FastFlowFlag=(Pos>0.5)+(Neg>0.5);
% FastPosFlow=FastFlowFlag.*(PosNegDiff<0.5).*(PosNegDiff>0);
% FastNegFlow=FastFlowFlag.*(PosNegDiff>-0.5).*(PosNegDiff<0);
% TrueFlow = FastFlowFlag+(abs(PosNegDiff)>0.25);
% Pos=(Pos+FastPosFlow.*Neg*2).*TrueFlow;
% Neg=(Neg+FastNegFlow.*Pos*2).*TrueFlow; 

%% Flow reconstruction by ratio between the linearly filtered image (FlowMap) and the original amplitude image
LIMPFlow = 1.414*pi*(Pos-Neg);

%% filtering in complex domain
% LIMPFlow = LIMPFlow .* (LIMPFlow <= pi) .* (LIMPFlow >= -pi) + pi*((LIMPFlow < -pi) + (LIMPFlow > pi));
% LIMPFlow = IntensityImage .* exp(1i*LIMPFlow);
% LIMPFlow = FiltImage(LIMPFlow, [3,3], 'Filter: mean');
% LIMPFlow = angle(LIMPFlow);

%% alternative LIMP1
% function [LIMPFlow,IntensityImage]=GenerateLIMPFlow(CaliedData, FullRange, LIMPFilterSize,LIMPFilterType, ZeroPadNum, MSS)
% 
% 
% [ZeroPadNum, AlineNum]=size(CaliedData);
% IntensityImage=FiltImage(abs(ifft(CaliedData,ZeroPadNum)),LIMPFilterSize,LIMPFilterType);
% if MSS
%     CaliedData=diff(CaliedData,1,2);
%     AlineNum=AlineNum-1;
%     IntensityImage=((IntensityImage(:,2:end))+(IntensityImage(:,1:end-1)))/2;
% end
% ComplexImage=ifft(CaliedData,ZeroPadNum);
% 
% % ComplexImage=FiltImage(ifft(CaliedData,ZeroPadNum), LIMPFilterSize,LIMPFilterType);
% % CaliedData=fft(ComplexImage);
% % CaliedData=CaliedData(1:ZeroPadNum,:);
% 
% if ~FullRange
%     ComplexImage=ComplexImage(2:end/2,:);
%     IntensityImage=IntensityImage(2:end/2,:);
% end
% 
% fftData=fft(CaliedData,[],2);
% 
% %% the kernel as Filter could be very important for quantify laminal flows
% Filter=zeros(1,AlineNum);
% Filter(2:end)=triang(AlineNum-1);
% Filter(end/2:end)=0;
% LIMPPosFlowMap=(ifft(ifft(fftData*sparse(diag(Filter)),[],2),ZeroPadNum));
% LIMPPosFlowMap=abs(FiltImage(LIMPPosFlowMap(2:end/2,:)./ComplexImage,LIMPFilterSize,LIMPFilterType));
% % LIMPPosFlowMap=abs(FiltImage(LIMPPosFlowMap,  );
% 
% Filter(2:end)=gausswin(AlineNum-1);
% Filter(2:end/2)=0;
% LIMPNegFlowMap=(ifft(ifft(fftData*sparse(diag(Filter)),[],2),ZeroPadNum));
% LIMPNegFlowMap=abs(FiltImage(LIMPNegFlowMap(2:end/2,:)./ComplexImage,LIMPFilterSize,LIMPFilterType));
% % LIMPNegFlowMap=abs(FiltImage(LIMPNegFlowMap,  LIMPFilterSize,LIMPFilterType));
% 
% PosMask=LIMPPosFlowMap>=LIMPNegFlowMap;%in order to resolve the part close to PI and 0
% NegMask=~PosMask;
% Filter=zeros(1,AlineNum);
% Filter(2:end)=triang(AlineNum-1);
% LIMPFlow=(ifft(ifft(fftData*sparse(diag(Filter)),[],2),ZeroPadNum));
% LIMPFlow=abs(FiltImage(LIMPFlow(2:end/2,:)./ComplexImage,LIMPFilterSize,LIMPFilterType));
% 
% LIMPFlow=pi*(LIMPFlow.*PosMask-LIMPFlow.*NegMask); %LIMPV=Diff
% LIMPFlow=pi*(LIMPPosFlowMap.*PosMask-LIMPNegFlowMap.*NegMask); %LIMPV=Diff

%% alternative LIMP2
% function [LIMPFlow,IntensityImage]=GenerateLIMPFlow(CaliedData, FullRange, LIMPFilterSize,LIMPFilterType, ZeroPadNum, MSS)
% 
% IntensityImage=FiltImage(abs(ifft(CaliedData,ZeroPadNum)),LIMPFilterSize,LIMPFilterType);
% if MSS
%     CaliedData=diff(CaliedData,1,2);
%     IntensityImage=((IntensityImage(:,2:end))+(IntensityImage(:,1:end-1)))/2;
% end
% 
% ComplexImage=ifft(CaliedData,ZeroPadNum);
% 
% if ~FullRange
%     ComplexImage=ComplexImage(2:end/2,:);
%     IntensityImage=IntensityImage(2:end/2,:);
% end
% 
% % IntensityImage=IntensityImage/median(max(IntensityImage))/10;
% 
% AlineNum=size(CaliedData,2);
% fftData=fft(CaliedData,[],2);
% 
% %% the kernel as Filter could be very important for quantify laminal flows
% Filter=zeros(1,AlineNum);
% Filter(2:AlineNum/2+1)=(1:AlineNum/2)/(AlineNum/2);
% 
% LIMPPosFlowMap=(ifft(ifft(fftData*sparse(diag(Filter)),[],2),ZeroPadNum));
% % LIMPPosFlowMap=abs(FiltImage(LIMPPosFlowMap,  LIMPFilterSize,LIMPFilterType))./abs(FiltImage(Intensity,  LIMPFilterSize,LIMPFilterType));
% LIMPPosFlowMap=LIMPPosFlowMap(2:end/2,:)./ComplexImage;
% LIMPPosFlowMap=abs(FiltImage(LIMPPosFlowMap,  LIMPFilterSize,LIMPFilterType));
% % LIMPPosFlowMap=min(LIMPPosFlowMap,1);
% % figure;imshow(LIMPPosFlowMap);colormap hot;
% 
% Filter=zeros(1,AlineNum);
% Filter(AlineNum:-1:AlineNum/2+1)=(1:AlineNum/2)/(AlineNum/2);
% LIMPNegFlowMap=(ifft(ifft(fftData*sparse(diag(Filter)),[],2),ZeroPadNum));
% % LIMPNegFlowMap=abs(LIMPNegFlowMap(2:end/2,:));
% % LIMPNegFlowMap=abs(FiltImage(LIMPNegFlowMap,  LIMPFilterSize,LIMPFilterType))./abs(FiltImage(Intensity,  LIMPFilterSize,LIMPFilterType)); 
% LIMPNegFlowMap=LIMPNegFlowMap(2:end/2,:)./ComplexImage;
% LIMPNegFlowMap=abs(FiltImage(LIMPNegFlowMap, LIMPFilterSize,LIMPFilterType)); 
% % LIMPNegFlowMap=min(LIMPNegFlowMap,1);   
% % clear fftData;
% % figure;imshow(LIMPNegFlowMap);colormap hot;
% 
% % PosMask=LIMPPosFlowMap>(LIMPNegFlowMap.*(LIMPRatio-(LIMPRatio-1)*LIMPNegFlowMap));%in order to resolve the part close to PI
% % NegMask=LIMPNegFlowMap>(LIMPPosFlowMap.*(LIMPRatio-(LIMPRatio-1)*LIMPPosFlowMap));
% 
% % PosMask=LIMPPosFlowMap>(LIMPNegFlowMap.*(Const('LIMPRatio')-(Const('LIMPRatio')-1)*LIMPNegFlowMap));%in order to resolve the part close to PI and 0
% % NegMask=LIMPNegFlowMap>(LIMPPosFlowMap.*(Const('LIMPRatio')-(Const('LIMPRatio')-1)*LIMPPosFlowMap));
% PosMask=LIMPPosFlowMap>=LIMPNegFlowMap;%in order to resolve the part close to PI and 0
% NegMask=~PosMask;

%% residual 
% % if VDiffOrMax
% %     LIMPFlowV=(LIMPPosFlowMap-LIMPNegFlowMap).*(PosMask+NegMask); %LIMPV=Diff
% % else
% %     LIMPFlowV=max(LIMPPosFlowMap,LIMPNegFlowMap).*(PosMask-NegMask);%LIMPV=larger one
% % end
% % 
% % if PSumOrMin
% %     LIMPFlowP=(LIMPPosFlowMap+LIMPNegFlowMap);%LIMPP=sum
% % else
% %     LIMPFlowP=min(LIMPPosFlowMap,LIMPNegFlowMap)*2; %LIMPP=smaller one in both side
% % end
% 
% %% calculate LIMPFlow as the sum of Vertical and Parallel flow
% % LIMPFlow=max(LIMPPosFlowMap,LIMPNegFlowMap).*(PosMask-NegMask)*pi; %LIMPV=Diff
% LIMPFlow=pi*(LIMPPosFlowMap.*PosMask-LIMPNegFlowMap.*NegMask); %LIMPV=Diff