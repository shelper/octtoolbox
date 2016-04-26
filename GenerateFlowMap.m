function FlowMap=GenerateFlowMap(CaliedData, ZeroPadNum, Window)
%generate angiographic flow imaging using DPRM

    AlineNum=size(CaliedData,2);
    ZeroPad=round(AlineNum*(1-diff(Window))/diff(Window));
    Filter1=gausswin(AlineNum/2,3.6)';Filter1(end/2:end)=1;Filter2=fliplr(Filter1);
    Filter=[Filter2,zeros(1,ZeroPad),Filter1];
    Filter=imresize(Filter, [1,AlineNum]);
    Filter=circshift(Filter,[1,round(mean(Window)*AlineNum)]);
    Filter=Filter+fliplr(Filter);
%     Filter=zeros(AlineNum,1);
%     Filter(Window(1)*end+1:Window(2)*end)=1;
    fftData=fft(CaliedData,[],2);
    CaliedData=ifft(fftData*sparse(diag(Filter)),[],2);
    FlowMap=abs(ifft(CaliedData,ZeroPadNum)*ZeroPadNum)+1;

    %% use weightedratio map and intensity map
%     PosHilbertImageRatio = PosHilbertImage(2:end/2,:)./PosHilbertImage(end:-1:end/2+2,:);
%     NegHilbertImageRatio = NegHilbertImage(end:-1:end/2+2,:)./NegHilbertImage(2:end/2,:);              
%     DFRMFlowMap=PosHilbertImageRatio+NegHilbertImageRatio; 
%     DFRMFlowMap=PosHilbertImage(2:end/2,:)+NegHilbertImage(end:-1:end/2+2,:); 
%     DFRMFlowMap=FiltImage(DFRMFlowMap,FilterSize,FilterType);

              