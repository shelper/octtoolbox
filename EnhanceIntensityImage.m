function EnhanceIntensityImage(handles,FilePath, FileName, ProcedureN)
global hIntensityImage IntensityImage EnhancedIntensityImage ScaledImage;

if ProcedureN<=1
    FreqFilterThresh=get(handles.FreqFilterThresh,'String');
    FreqFilterThresh=str2double(FreqFilterThresh(12:end));    
    FreqIntensityThresh=get(handles.FreqIntensityThresh,'String');
    FreqIntensityThresh=str2double(FreqIntensityThresh(12:end)); 
    EnhancedIntensityImage=LateralFreqFilter(IntensityImage,FreqFilterThresh,FreqIntensityThresh);

    FilterType=get(handles.FilterType,'String');
    FilterType=FilterType{get(handles.FilterType,'Value')};
    FilterSize=get(handles.FilterSize,'String');
    FilterSize=str2num(FilterSize(13:end));
    EnhancedIntensityImage=FiltImage(EnhancedIntensityImage, FilterSize, FilterType);

    ResizeRatio=get(handles.ImgResize,'String');
    ResizeRatio=str2num(ResizeRatio(9:end));
    NewSize=ResizeRatio.*size(EnhancedIntensityImage);
    EnhancedIntensityImage=imresize(EnhancedIntensityImage,NewSize,'box');
end

if ProcedureN<=2
    ImageScale=get(handles.ImageScale,'String');
    switch ImageScale{get(handles.ImageScale,'Value')};
        case 'Log Scale'
            ScaledImage=20*log10(EnhancedIntensityImage);
        case 'Linear Scale'
            ScaledImage=EnhancedIntensityImage;      
    end
end

ImageFOV=get(handles.ImageFOV,'String');
ImageFOV=str2num(ImageFOV(5:end));          
DisplayedIntensityImage=ScaledImage(ImageFOV(1)*end+1:ImageFOV(2)*end,ImageFOV(3)*end+1:ImageFOV(4)*end);
if get(handles.BackForward,'Value')
    DisplayedIntensityImage=DisplayedIntensityImage(:,1:end/2)+DisplayedIntensityImage(:,end/2+1:end);
end    

DisplayRange=get(handles.DisplayRange,'String');
DisplayRange=GetDisplayRange(DisplayedIntensityImage, str2num(DisplayRange(7:end)));
DisplayedIntensityImage=(DisplayedIntensityImage-DisplayRange(1))/(DisplayRange(2)-DisplayRange(1));

if ishandle(hIntensityImage)
    set(hIntensityImage,'CData', DisplayedIntensityImage);       
else
    figure('Name','IntensityImage','NumberTitle','off');
    hIntensityImage=imshow(DisplayedIntensityImage);colormap gray;impixelinfo;
end        

%% save
if get(handles.AutoSaveIntensityImage,'Value')
%     fid = fopen([FilePath,FileName(4:end)],'w');fwrite(fid, EnhancedIntensityImage, 'float');fclose(fid);
    imwrite(uint16(DisplayedIntensityImage*65535),[FilePath,FileName(1:end-4),'Flt.tif']);
end
%% save the frame diff image for further angiographic study
% IterationNum=8;
% SurfaceDepth=80;
% PenatratedDepth=500;
% Thickness=PenatratedDepth-SurfaceDepth+1;
% IntensityImage=IntensityImage(SurfaceDepth:PenatratedDepth,:);%this one has problem because only one time cropping is needed.
% ScanNum=size(IntensityImage,2)/IterationNum;
% if get(handles.AngioByStd,'Value') && (FileName(end-4)=='B')
%     StackImage=zeros(Thickness,ScanNum-ScanNum/10,IterationNum);
%     TempImage=abs(EnhancedIntensityImage(:,ScanNum/10+1:ScanNum));%pls input the
%     StackImage(:,:,1)=TempImage;
%     for i=1:IterationNum-1
%         SubImage=EnhancedIntensityImage(:,ScanNum*(i+0.1)+1:ScanNum*(i+1));
%         ImageCC=normxcorr4smallshift(TempImage,SubImage,5, 5);
%         [MaxCC,Index]=max(ImageCC(:));
%         [RowShift,ColShift]=ind2sub(size(ImageCC),Index(1));
%         RowShift=RowShift-6;
%         ColShift=ColShift-6;
%         StackImage(:,:,i+1)=circshift(SubImage,[RowShift,ColShift]);
%         StackImage(:,:,i+1)=SubImage;
%     end
%     EnhancedIntensityImage=mean(StackImage,3);
%     
%     FlowMap=std(StackImage,1,3)./EnhancedIntensityImage.*(EnhancedIntensityImage>0.25);
%     if ishandle(hFlowMapImage)
%         set(hFlowMapImage,'CData', FlowMap);       
%     else
%         figure('Name','FlowMap','NumberTitle','off');
%         hFlowMapImage=imshow(FlowMap);colormap hot;impixelinfo;
%     end    
%     imwrite(uint16(FlowMap*65535),[FilePath,FileName(1:end-4),'StdMap','.tif']); 
% 
% elseif (FileName(end-4)=='F')
%         EnhancedIntensityImage=imresize(EnhancedIntensityImage,size(EnhancedIntensityImage)./[1,IterationNum],'box');
%         EnhancedIntensityImage=EnhancedIntensityImage(:,ScanNum/10+1:ScanNum);
% end
