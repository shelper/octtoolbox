function ReconstructDFRMFlow(handles,ZeroPadNum, FilePath, FileName, ProcedureN)
global CaliedData hDFRMFlowImage IntensityImage DFRMFlow ThreshMask EdgeMask DemotionedDFRMFlow EnhancedDFRMFlow;

if ProcedureN==0
    %% reconstruction
    DFRMFlowResolution=get(handles.DFRMFlowResolution,'String');
    DFRMFlowResolution=str2num(DFRMFlowResolution(13:end));

    DFRMFilterType=get(handles.DFRMFilterType,'String');
    DFRMFilterType=DFRMFilterType{get(handles.DFRMFilterType,'Value')};%check if the filter name is correct;
    DFRMFilterSize=get(handles.DFRMFilterSize,'String');
    DFRMFilterSize=str2num(DFRMFilterSize(13:end));

    GPU=get(handles.GPU,'String');
    GPU=GPU{get(handles.GPU,'Value')};
    GPU=strcmp('Run on GPU',GPU);
    DFRMFlow = GenerateDFRMFlow(CaliedData, ZeroPadNum, DFRMFlowResolution, DFRMFilterType, DFRMFilterSize, GPU);

end
%% Image Enhancing
if ProcedureN<=1
    DFRMFlow=DFRMFlow-2*(DFRMFlow>1);

    FilterType=get(handles.FilterType,'String');
    FilterType=FilterType{get(handles.FilterType,'Value')};
    FilterSize=get(handles.FilterSize,'String');
    FilterSize=str2num(FilterSize(13:end));
    EnhancedDFRMFlow=FiltFlowImage(DFRMFlow, FilterSize, FilterType);

end

if ProcedureN<=2
    IntensityThresh=get(handles.ThresholdDFRM,'String');
    IntensityThresh=str2double(IntensityThresh(12:end));
    if IntensityThresh
        ThreshMask=im2bw(IntensityImage,IntensityThresh);
    else
        ThreshMask=ones(size(DFRMFlow));
    end
    ThreshMask(1:end/20,:)=0;
end


%% demotion
if ProcedureN<=3
    ImageFOV=get(handles.ImageFOV,'String');
    ImageFOV=str2num(ImageFOV(5:end));          
    DemotionedDFRMFlow=(EnhancedDFRMFlow(ImageFOV(1)*end+1:ImageFOV(2)*end,ImageFOV(3)*end+1:ImageFOV(4)*end)...
        .*ThreshMask(ImageFOV(1)*end+1:ImageFOV(2)*end,ImageFOV(3)*end+1:ImageFOV(4)*end));
    if get(handles.BackForward,'Value')
        tempDivArray=ones(size(DemotionedDFRMFlow).*[1,0.5])+(DemotionedDFRMFlow(:,1:end/2)~=0).*(DemotionedDFRMFlow(:,end/2+1:end)~=0);
        DemotionedDFRMFlow=(DemotionedDFRMFlow(:,1:end/2)-DemotionedDFRMFlow(:,end/2+1:end))./tempDivArray;clear tempArray;    
    end    
    if get(handles.DFRMDemotion,'Value')
        DemotionedDFRMFlow= RemoveBulkyPhaseNoise(DemotionedDFRMFlow);
    end
end

ResizeRatio=get(handles.ImgResize,'String');
ResizeRatio=str2num(ResizeRatio(9:end));
NewSize=ResizeRatio.*size(DemotionedDFRMFlow);
DisplayedDFRMFlow=imresize(DemotionedDFRMFlow,NewSize,'box');
DisplayedDFRMFlow=uint16(DisplayedDFRMFlow*32768+32768);
if ishandle(hDFRMFlowImage)
    set(hDFRMFlowImage,'CData',DisplayedDFRMFlow(100:600,:));       
else
    figure('Name','DFRMFlow','NumberTitle','off');
    hDFRMFlowImage=imshow(DisplayedDFRMFlow);colormap jet;impixelinfo;
end
%% save
if get(handles.AutoSaveDFRMFlow,'Value')
    imwrite(DisplayedDFRMFlow,[FilePath,FileName(1:end-4),'DFRM.tif']);
end
% end
   
