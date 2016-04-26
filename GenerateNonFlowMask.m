

function NonFlowMask=GenerateNonFlowMask(Flow,SpanSize,VarThresh, DilateSize, CloseSize, ErodeSize)

%% correct the zero points caused var change
SpanWindow=ones(SpanSize,1)/SpanSize;
CorrGain=imfilter(double(Flow~=0),SpanWindow);
Flow=Flow+(Flow==0).*imfilter(Flow,SpanWindow)./CorrGain;

VarImg=imfilter(Flow.^2,SpanWindow)-(imfilter(Flow,SpanWindow).^2);
VarImg=VarImg./CorrGain;

se=strel('line',DilateSize,0);
NonFlowMask=imdilate(VarImg>VarThresh,se);

se=strel('disk',CloseSize);
NonFlowMask=imclose(NonFlowMask,se);

se=strel('disk',ErodeSize);
NonFlowMask=imerode(NonFlowMask,se);

NonFlowMask=(~NonFlowMask);