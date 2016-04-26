function DFRMFlow=GenerateDFRMFlow(CaliedData, ZeroPadNum, DFRMFlowResolution, DFRMFilterType, DFRMFilterSize, GPU)
% global FilePath FileName;
DFRMProgressBar = waitbar(0,'DFRM Quantitative measurement...');
%% USE GPU
if GPU 
    PosBinaryFlowRatio=gsingle(0);
    NegBinaryFlowRatio=gsingle(0);
    for PhaseShift=-0.5:DFRMFlowResolution:0.5
        RampedData = RampPhase(CaliedData,PhaseShift); 
        HilbertImage=abs(ifft(hilbert(gsingle(RampedData'))',ZeroPadNum));
        % HilbertImage=FiltImage(HilbertImage,gsingle(DFRMFilterSize),DFRMFilterType);
        % Calculate the Postive Flow using piling up of the binarized existance
        if PhaseShift<0
            NegBinaryFlowRatio=NegBinaryFlowRatio+(HilbertImage(end:-1:end/2+2,:)>HilbertImage(2:end/2,:)); 
        else
            PosBinaryFlowRatio=PosBinaryFlowRatio+(HilbertImage(2:end/2,:)>HilbertImage(end:-1:end/2+2,:)); 
        end
%         clear HilbertImage;
        waitbar(PhaseShift+0.5,DFRMProgressBar);    
    end
    DFRMFlow = DFRMFlowResolution*double(((PosBinaryFlowRatio>NegBinaryFlowRatio).*(PosBinaryFlowRatio+NegBinaryFlowRatio-0.5) ...
        -(NegBinaryFlowRatio>=PosBinaryFlowRatio).*(PosBinaryFlowRatio+NegBinaryFlowRatio-0.5)));
%% USE CPU
else
    PosBinaryFlowRatio=0;
    NegBinaryFlowRatio=0;
%     i=1;
    for PhaseShift=-0.5:DFRMFlowResolution:0.5
        RampedData = RampPhase(CaliedData,PhaseShift); 
        HilbertImage=abs(ifft(hilbert(RampedData')', ZeroPadNum));
        HilbertImage=FiltImage(HilbertImage,DFRMFilterSize,DFRMFilterType);        
        % Calculate the Postive Flow using piling up of the binarized existance
        if PhaseShift<0
            NegBinaryFlowRatio=NegBinaryFlowRatio+(HilbertImage(end:-1:end/2+2,:)>HilbertImage(2:end/2,:)); 
        else
            PosBinaryFlowRatio=PosBinaryFlowRatio+(HilbertImage(2:end/2,:)>HilbertImage(end:-1:end/2+2,:)); 
        end
%         %% to save intensity image from each ramping step  
%             HilbertImage=HilbertImage(end/100:end/2,:);
%             HilbertImage=uint16(HilbertImage/max(HilbertImage(:))*65535);
%             imwrite(HilbertImage(end/100:end/2,:),[FilePath,FileName(1:end-4),num2str(PhaseShift),'.tif'],'Compression','none');
%             Tmp(i)=HilbertImage(177,580);
%             i=i+1;

        waitbar(PhaseShift+0.5,DFRMProgressBar);    
    end
    

    DFRMFlow = DFRMFlowResolution*((PosBinaryFlowRatio>NegBinaryFlowRatio).*(PosBinaryFlowRatio+NegBinaryFlowRatio-0.5) ...
        -(NegBinaryFlowRatio>=PosBinaryFlowRatio).*(PosBinaryFlowRatio+NegBinaryFlowRatio-0.5));
end

DFRMFlow=DFRMFlow-2*(DFRMFlow>1);
close(DFRMProgressBar);
