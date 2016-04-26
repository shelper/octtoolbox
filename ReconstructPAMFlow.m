function CaliedData=ReconstructPAMFlow(CaliedData,PAMPar,ImgPar,FlowPar,FileName, ProcedureN)
global hPAMFlowImage PAMFlow DemotionedPAMFlow IntensityMask ;

% use procedure number ProcedureN to identify which controller is triggered, thus save processing time.
if ProcedureN==0
    %% reconstruction
    ComplexImage=ifft(CaliedData,ImgPar.ZeroPadNum);
%     ComplexImage=ComplexImage(2:end/2,:);    
    PAMFlow = zeros(size(ComplexImage));
    PAMFlow(:,2:end) = angle(ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1)))/pi;
    PAMFlow(:,1) = angle(ComplexImage(:,1))/pi; 
    if ImgPar.fftShift
        PAMFlow=fftshift(PAMFlow,1);
    end
%     PAMFlow(1:50,:) =0;
    
end


%%  FOV & Thresh & demotion%%
if ProcedureN<=1
    DemotionedPAMFlow=PAMFlow(ImgPar.ImageFOV(1)*end+1:ImgPar.ImageFOV(2)*end,ImgPar.ImageFOV(3)*end+1:ImgPar.ImageFOV(4)*end);
    DemotionedPAMFlow=DemotionedPAMFlow.*IntensityMask;
    if FlowPar.RemoveBMA
        NonFlowMask=GenerateNonFlowMask(DemotionedPAMFlow, ...
            FlowPar.BMAStdSize,FlowPar.BMAStdThresh, FlowPar.BMADilate, FlowPar.BMAClose, FlowPar.BMAErode);
        NonFlowMask=1;
        PAMBMA= GenBulkyPhaseNoise(DemotionedPAMFlow.*NonFlowMask);
%         figure;imshow(NonFlowMask);
%         clear NonFlowMask;
        DemotionedPAMFlow = angle(exp(1i*DemotionedPAMFlow*pi)*sparse(diag(exp(-1i*PAMBMA*pi)))).*(DemotionedPAMFlow~=0)/pi;
    end
 end
%% Filtering and resize
if ProcedureN<=2
    EnhancedPAMFlow=FiltFlowImage(DemotionedPAMFlow, ImgPar.ImgFilterSize, ImgPar.ImgFilterType);    
    EnhancedPAMFlow=imresize(EnhancedPAMFlow,ImgPar.ImgResize.*size(EnhancedPAMFlow),'box');
    if FlowPar.AbsFlow
        EnhancedPAMFlow=abs(EnhancedPAMFlow);
        EnhancedPAMFlow=uint16(EnhancedPAMFlow*65535);        
    else
        EnhancedPAMFlow=uint16(EnhancedPAMFlow*32768+32768);
    end
end

%% Display and save
if ProcedureN<=3
    if ishandle(hPAMFlowImage)
        set(hPAMFlowImage,'CData',EnhancedPAMFlow);       
    else
        figure('Name','PAMFlow','NumberTitle','off');
        hPAMFlowImage=imshow(EnhancedPAMFlow);colormap jet;impixelinfo;
    end
    %% save
    if PAMPar.AutoSave
        imwrite(EnhancedPAMFlow,[FileName(1:end-4),'PAM.tif']);
    end
end

%% Generate Demotioned CaliedData
if nargout==1 
        PAMBMASum=PAMBMA*triu(ones(size(PAMBMA,2)));
        CaliedData(:,ImgPar.ImageFOV(3)*end+1:ImgPar.ImageFOV(4)*end)=CaliedData(:,ImgPar.ImageFOV(3)*end+1:ImgPar.ImageFOV(4)*end)*sparse(diag(exp(-1i*PAMBMASum*pi)));
end


