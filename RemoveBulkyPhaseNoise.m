function [DenoisedFlow,BMA]=RemoveBulkyPhaseNoise(Flow)
AlineFlowStd1=nonzerostd(Flow);
AlineFlowStd2=nonzerostd(Flow+2*(Flow<0));
WrapLines = (AlineFlowStd1>AlineFlowStd2);
Flow=Flow*sparse(diag(~WrapLines))+(Flow+(Flow<0)*2)*sparse(diag(WrapLines));

tmpFlow = sort(Flow+(Flow==0)*4);% let flow==0 to be excluded in sorting by setting them to be largest
BMA = tmpFlow( round(sum(Flow~=0)/2) + (1:size(Flow,1):numel(Flow)) );
DenoisedFlow = angle(exp(1i*Flow*pi)*sparse(diag(exp(-1i*BMA*pi)))).*(Flow~=0)/pi;
end

function [DenoisedFlow,BMA]=RemoveBulkyPhaseNoise6(Flow)
AlineFlowStd1=nonzerostd(Flow);
AlineFlowStd2=nonzerostd(Flow+2*(Flow<0));
WrapLines = (AlineFlowStd1>AlineFlowStd2);
Flow=Flow*sparse(diag(~WrapLines))+(Flow+(Flow<0)*2)*sparse(diag(WrapLines));

tmpFlow = sort(Flow+(Flow==0)*4);% let flow==0 to be excluded in sorting by setting them to be largest
BMA = tmpFlow( round(sum(Flow~=0)/2) + (1:size(Flow,1):numel(Flow)) );
DenoisedFlow = angle(exp(1i*Flow*pi)*sparse(diag(exp(-1i*BMA*pi)))).*(Flow~=0)/pi;
end


function [DenoisedFlow,BMA]=RemoveBulkyPhaseNoise5(Flow)
% StdThresh=0.1;
% 
% AlineNum=size(Flow,2);
% % AlineFlowStd1=nonzerostd(Flow);
% % AlineFlowStd2=nonzerostd(Flow+2*(Flow<0));
% % WrapLines = (AlineFlowStd1>AlineFlowStd2);
% % Flow=Flow*sparse(diag(~WrapLines))+(Flow+(Flow<0)*2)*sparse(diag(WrapLines));
% 
% BMA=zeros(1,AlineNum);
% [foo,AlineIndex]=min(nonzerostd(Flow));
% 
% Aline=Flow(:,AlineIndex);
% BMA(AlineIndex)=median(Aline(Aline~=0));
% AlineMask=(Aline<BMA(AlineIndex)+StdThresh).*(Aline>BMA(AlineIndex)-StdThresh);
% for n=AlineIndex-1:-1: 1
%     Aline=Flow(:,n);
%     AlineMask=AlineMask.*(Aline~=0);
%     se = strel('disk',3);
%     AlineMask = imopen(AlineMask,se);
%     AlineMask=AlineMask+(sum(AlineMask)==0);
%     Aline=Aline+2*(std(Aline(logical(AlineMask)))>std(Aline(logical(AlineMask))+2));
%     BMA(n)=median(Aline(logical(AlineMask)));    
% %     if n==410
% %         figure;plot(AlineMask);
% %         hold on; plot(Aline)
% %     end
% %     
%     AlineMask = (Aline<BMA(n)+StdThresh).*(Aline>BMA(n)-StdThresh);
% end
% 
% Aline=Flow(:,AlineIndex);
% AlineMask = (Aline<BMA(AlineIndex)+StdThresh).*(Aline>BMA(AlineIndex)-StdThresh);
% for n=AlineIndex+1:AlineNum 
%     Aline=Flow(:,n);
%     AlineMask=AlineMask.*(Aline~=0);
%     se = strel('disk',3);
%     AlineMask = imopen(AlineMask,se);
%     AlineMask=AlineMask+(sum(AlineMask)==0);
%     Aline=Aline+2*(std(Aline(logical(AlineMask)))>std(Aline(logical(AlineMask))+2));
%     BMA(n)=median(Aline(logical(AlineMask)));    
% %     if n==410
% %         figure;plot(AlineMask);
% %         hold on; plot(Aline)
% %     end
% 
%     AlineMask = (Aline<BMA(n)+StdThresh).*(Aline>BMA(n)-StdThresh);
% end
% 
% DenoisedFlow = angle(exp(1i*Flow*pi)*sparse(diag(exp(-1i*BMA*pi)))).*(Flow~=0)/pi;
end

function [DenoisedFlow,BMA]=RemoveBulkyPhaseNoise4(Flow)
% AlineFlowStd1=nonzerostd(Flow);
% AlineFlowStd2=nonzerostd(Flow+2*(Flow<0));
% WrapLines = (AlineFlowStd1>AlineFlowStd2);
% Flow=Flow*sparse(diag(~WrapLines))+(Flow+(Flow<0)*2)*sparse(diag(WrapLines));
% 
% tmpFlow = sort(Flow+(Flow==0)*4);% let flow==0 to be excluded in sorting by setting them to be largest
% BMA = tmpFlow( round(sum(Flow~=0)/2) + (1:size(Flow,1):numel(Flow)) );
% DenoisedFlow = angle(exp(1i*Flow*pi)*sparse(diag(exp(-1i*BMA*pi)))).*(Flow~=0)/pi;
end

function DenoisedFlow=RemoveBulkyPhaseNoise1(Flow)
% for N=1:2 %iternation number2
%     DenoisedComplexImage=ComplexImage./AvgComplexImage;
%     NegAngleComplexImage=ComplexImage.*(angle(DenoisedComplexImage)>0.1);
%     PosAngleComplexImage=ComplexImage.*(angle(DenoisedComplexImage)<-0.1);
% 
%     %select the A-lines that has no major blood vessels
% %     NumNeg=sum(NegAngleComplexImage~=0);
% %     NumPos=sum(PosAngleComplexImage~=0);
% %     NumTotal=sum(FlowMask);
% %     BigFlowLine=((NumNeg./NumTotal)>0.9) .* ((NumPos./NumTotal)>0.9);
% %     ComplexImageNoBigVessel=ComplexImage.*repmat(~BigFlowLine,[size(ComplexImage,1),1]);
% 
%     %remove the flow part of the A-line that has major blood vessels
% %     NegLine=(nonzerostd(NegAngleComplexImage)<nonzerostd(PosAngleComplexImage)).*BigFlowLine;
% %     PosLine=BigFlowLine-NegLine;
% %     ComplexImageNoBigVessel=ComplexImageNoBigVessel+NegAngleComplexImage.*repmat(NegLine,[size(ComplexImage,1),1]) ...
% %                                                    +PosAngleComplexImage.*repmat(PosLine,[size(ComplexImage,1),1]);
%     PosLine=((nonzerostd(NegAngleComplexImage)-nonzerostd(PosAngleComplexImage))>0.3);
%     NegLine=((nonzerostd(PosAngleComplexImage)-nonzerostd(NegAngleComplexImage))>0.3);
%     OtherLine=1-PosLine-NegLine;
%     ComplexImageNoBigVessel=NegAngleComplexImage.*repmat(NegLine,[size(ComplexImage,1),1]) ...
%         +PosAngleComplexImage.*repmat(PosLine,[size(ComplexImage,1),1]) ...
%         +ComplexImage.*repmat(OtherLine,[size(ComplexImage,1),1]) ;
% 
%     AvgComplexImage=repmat(sum(ComplexImageNoBigVessel),[size(ComplexImageNoBigVessel,1),1]);
%     AvgComplexImage=AvgComplexImage+(AvgComplexImage==0);
% end
end

function PhaseNoise=RemoveBulkyPhaseNoise2(Flow,FlowMask)
% FilteredFlow=FiltImage(Flow,[10,1],'Filter: mean',1);
% % FilteredFlow=Flow;
% ComplexImage=cos(FilteredFlow*pi)+ (1i*sin(FilteredFlow*pi));
% ComplexImage=ComplexImage.*FlowMask;
% 
% AvgComplexImage=repmat(sum(ComplexImage),[size(ComplexImage,1),1]);
% DenoisedComplexImage=ComplexImage./AvgComplexImage;
% NegAngleComplexImage=ComplexImage.*(angle(DenoisedComplexImage./AvgComplexImage)>0.1);
% PosAngleComplexImage=ComplexImage.*(angle(DenoisedComplexImage./AvgComplexImage)<0.1);
% 
% NumNeg=sum(NegAngleComplexImage~=0);
% NumPos=sum(PosAngleComplexImage~=0);
% NumTotal=sum(FlowMask);
% BigFlowLine=((NumNeg./NumTotal)>0.2) .* ((NumPos./NumTotal)>0.2);
% FilteredComplexImage=ComplexImage.*repmat(~BigFlowLine,[size(ComplexImage,1),1]);
% 
% NegLine=(nonzerostd(NegAngleComplexImage)<nonzerostd(PosAngleComplexImage)).*BigFlowLine;
% PosLine=BigFlowLine-NegLine;
% FilteredComplexImage=FilteredComplexImage+NegAngleComplexImage.*repmat(NegLine,[size(ComplexImage,1),1]);
% FilteredComplexImage=FilteredComplexImage+PosAngleComplexImage.*repmat(PosLine,[size(ComplexImage,1),1]);
% 
% PhaseNoise=angle(sum(FilteredComplexImage))/pi;
end

function DenoisedFlow=RemoveBulkyPhaseNoise3(Flow,FlowMask)

% % FilteredFlow=FiltFlowImage(Flow,[1,1],'Filter: mean');
% % % FilteredFlow=Flow;
% % FilteredFlow(300:end,:)=0;
% ComplexImage=cos(FilteredFlow*pi)+ (1i*sin(FilteredFlow*pi));
% ComplexImage=ComplexImage.*FlowMask;
% 
% AvgComplexImage=repmat(sum(ComplexImage),[size(ComplexImage,1),1]);
% AvgComplexImage=AvgComplexImage+(AvgComplexImage==0);
% ComplexImage=(cos(Flow*pi)+ (1i*sin(Flow*pi)));
% DenoisedFlow=angle(ComplexImage./AvgComplexImage)/pi;

end

function MotionPhaseNoise=CaculateMotionPhaseNoise(Flow,FlowMask,NoiseThreshold)
% Flow=FiltFlowImage(Flow,[10,1],'mean');%TBF-? why [10,1]?, and shoudl filt before or after flow EdgeMask?
% SinFlow=sin(Flow*pi);
% CosFlow=cos(Flow*pi);
% Flow(:)=1;
% PointsNum=sum(FlowMask);
% PointsNum=PointsNum+(PointsNum==0);
% SinNoise2=sum(SinFlow.*FlowMask)./PointsNum;
% CosNoise2=sum(CosFlow.*FlowMask)./PointsNum;
% SinNoise1=0;
% CosNoise1=0;
% MotionPhaseNoise1=angle(CosNoise2+i*SinNoise2)/pi;
% 
% % figure;plot(angle(CosNoise2+i*SinNoise2)/pi,'r');
% 
% while sum((SinNoise2-SinNoise1).^2+(CosNoise2-CosNoise1).^2) > NoiseThreshold
%     SinNoise1=SinNoise2;
%     CosNoise1=CosNoise2;
%     %make 2 (low and high) 2D noise matrix to threshold the Background
%     SinFlowNoise=Flow*sparse(diag(SinNoise2)).*FlowMask;
%     CosFlowNoise=Flow*sparse(diag(CosNoise2)).*FlowMask;
%     NoiseMask= ((SinFlow-SinFlowNoise).^2+(CosFlow-CosFlowNoise).^2) < NoiseThreshold ;    
% 
%     PointsNum=sum(NoiseMask);
%     PointsNum=PointsNum+(PointsNum==0);
%     SinNoise2=sum(SinFlow.*NoiseMask)./PointsNum;
%     CosNoise2=sum(CosFlow.*NoiseMask)./PointsNum;
% %     hold on;plot(angle(CosNoise2+i*SinNoise2)/pi);
% end
% 
% MotionPhaseNoise2=angle(CosNoise2+i*SinNoise2)/pi;
% MotionPhaseNoise=MotionPhaseNoise2+MotionPhaseNoise1.*(MotionPhaseNoise2==0);

end
