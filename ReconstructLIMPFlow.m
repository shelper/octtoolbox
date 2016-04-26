function ReconstructLIMPFlow(CaliedData,LIMPPar,ImgPar,FlowPar,FileName, ProcedureN)
global hLIMPFlowImage LIMPFlow DemotionedLIMPFlow IntensityMask;

if ProcedureN==0
    LIMPFlow =GenerateLIMPFlow(CaliedData, ImgPar.ZeroPadNum, ...
        LIMPPar.LIMPFilterSize,LIMPPar.LIMPFilterType,LIMPPar.LIMPRatio, LIMPPar.PLIMPThresh,LIMPPar.VDiffOrMax, LIMPPar.PSumOrMin); 
    if ImgPar.fftShift
        LIMPFlow=fftshift(LIMPFlow,1);
    end
end

%%  FOV & Thresh & demotion%%
if ProcedureN<=1
    DemotionedLIMPFlow=LIMPFlow(ImgPar.ImageFOV(1)*end+1:ImgPar.ImageFOV(2)*end,ImgPar.ImageFOV(3)*end+1:ImgPar.ImageFOV(4)*end);
    DemotionedLIMPFlow=DemotionedLIMPFlow.*IntensityMask;
    if LIMPPar.LIMPDemotion
        NonFlowMask=GenerateNonFlowMask(DemotionedLIMPFlow, ...
            FlowPar.BMAStdSize,FlowPar.BMAStdThresh, FlowPar.BMADilate, FlowPar.BMAClose, FlowPar.BMAErode);
        LIMPBMA= GenBulkyPhaseNoise(DemotionedLIMPFlow.*NonFlowMask);
%         figure;imshow(NonFlowMask);
        DemotionedLIMPFlow = angle(exp(1i*DemotionedLIMPFlow*pi)*sparse(diag(exp(-1i*LIMPBMA*pi)))).*(DemotionedLIMPFlow~=0)/pi;
    end
end

%% Filtering and resize
if ProcedureN<=2
    EnhancedLIMPFlow=FiltFlowImage(DemotionedLIMPFlow, ImgPar.ImgFilterSize, ImgPar.ImgFilterType);    
%         EnhancedLIMPFlow=FiltImage(EnhancedLIMPFlow, ImgPar.ImgFilterSize, ImgPar.ImgFilterType);    
    EnhancedLIMPFlow=imresize(EnhancedLIMPFlow,ImgPar.ImgResize.*size(EnhancedLIMPFlow),'box');
    if FlowPar.AbsFlow
        EnhancedLIMPFlow=abs(EnhancedLIMPFlow);
        EnhancedLIMPFlow=uint16(EnhancedLIMPFlow*65535);        
    else
        EnhancedLIMPFlow=uint16(EnhancedLIMPFlow*32769+32768);
    end
end

%% Display and save
if ProcedureN<=3
    if ishandle(hLIMPFlowImage)
        set(hLIMPFlowImage,'CData',EnhancedLIMPFlow);       
    else
        figure('Name','LIMPFlow','NumberTitle','off');
        hLIMPFlowImage=imshow(EnhancedLIMPFlow);colormap jet;impixelinfo;
    end
    %% save
    if LIMPPar.AutoSave
        imwrite(EnhancedLIMPFlow,[FileName(1:end-4),'LIMP.tif']);
    end
end




   
