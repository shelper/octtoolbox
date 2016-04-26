function ReconstructDFRMFlowMap(handles, ZeroPadNum, FilePath, FileName, ProcedureN)
global CaliedData hDFRMFlowMapImage IntensityImage DFRMFlowMap ThreshMask EdgeMask EnhancedDFRMFlowMap;

% use procedure number ProcedureN to identify which controller is
% triggered, thus save processing time.
if ProcedureN==0
    %% reconstruction
    DFRMFlowMapSensitivity=get(handles.DFRMFlowMapSensitivity,'String');
    DFRMFlowMapSensitivity=str2double(DFRMFlowMapSensitivity(13:end));
    DFRMFlowMap = GenerateDFRMFlowMap(CaliedData, ZeroPadNum , DFRMFlowMapSensitivity);
end%% Image Enhancing

if ProcedureN<=1
    FilterType=get(handles.FilterType,'String');
    FilterType=FilterType{get(handles.FilterType,'Value')};
    FilterSize=get(handles.FilterSize,'String');
    FilterSize=str2num(FilterSize(13:end));
    EnhancedDFRMFlowMap=FiltImage(DFRMFlowMap, FilterSize, FilterType);
end


%% threshold
if ProcedureN<=2
    IntensityThresh=get(handles.ThresholdDFRMMap,'String');
    IntensityThresh=str2double(IntensityThresh(12:end));
    if IntensityThresh
        ThreshMask=im2bw(IntensityImage,IntensityThresh);
    else
        ThreshMask=ones(size(DFRMFlowMap));
    end
    ThreshMask(1:end/20,:)=0;
end

ImageFOV=get(handles.ImageFOV,'String');
ImageFOV=str2num(ImageFOV(5:end));          
DisplayedDFRMFlowMap=EnhancedDFRMFlowMap(ImageFOV(1)*end+1:ImageFOV(2)*end,ImageFOV(3)*end+1:ImageFOV(4)*end) ...
            .*ThreshMask(ImageFOV(1)*end+1:ImageFOV(2)*end,ImageFOV(3)*end+1:ImageFOV(4)*end);

if get(handles.BackForward,'Value')
    tempDivArray=ones(size(DisplayedDFRMFlowMap).*[1,0.5])+(DisplayedDFRMFlowMap(:,1:end/2)~=0).*(DisplayedDFRMFlowMap(:,end/2+1:end)~=0);
    DisplayedDFRMFlowMap=(DisplayedDFRMFlowMap(:,1:end/2)-DisplayedDFRMFlowMap(:,end/2+1:end))./tempDivArray;clear tempArray;    
end    

ResizeRatio=get(handles.ImgResize,'String');
ResizeRatio=str2num(ResizeRatio(9:end));
NewSize=ResizeRatio.*size(DisplayedDFRMFlowMap);

DisplayedDFRMFlowMap=imresize(DisplayedDFRMFlowMap,NewSize,'box');
DFRMFlowMapGain=get(handles.DFRMFlowMapGain,'String');
DFRMFlowMapGain=str2double(DFRMFlowMapGain(6:end));
DisplayedDFRMFlowMap=uint8(DisplayedDFRMFlowMap*DFRMFlowMapGain*255);
ColorMap=hot(256);
if ishandle(hDFRMFlowMapImage)
    set(hDFRMFlowMapImage,'CData',DisplayedDFRMFlowMap);       
else
    figure('Name','DFRMFlowMap','NumberTitle','off');
    hDFRMFlowMapImage=imshow(DisplayedDFRMFlowMap,ColorMap);impixelinfo;
end

%% save
if get(handles.AutoSaveDFRMFlowMap,'Value')
    imwrite(DisplayedDFRMFlowMap,ColorMap,[FilePath,FileName(1:end-4),'DFRMMap.tif']);
end
   


