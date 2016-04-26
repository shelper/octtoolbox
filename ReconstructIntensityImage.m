function ReconstructIntensityImage(CaliedData,PrePar, IntPar,ImgPar,FilePath, FileName, Procedure)
global hIntensityImage IntensityImage EnhancedIntensityImage AlineProfile;
debug=0;
if Procedure<=Const('Recon');
    IntensityImage=abs(fft(CaliedData,ImgPar.ZeroPadNum))+1;
    if strcmp(IntPar.ImageScale,'Log Scale');
        IntensityImage=20*log10(IntensityImage+1);
%         IntensityImage=IntensityImage-);
    end 
    
    if ~strcmpi(PrePar.DeConjugate,'None');
        IntensityImage=fftshift(IntensityImage,1);
    else
        IntensityImage=IntensityImage(1: end/2+2,:);
    end
end

IntensityImage=RemoveFreq(IntensityImage', 82:86, 2)';

if debug
%     CmpImg1=fft(CaliedData);
%     CmpImg1=fftshift(CmpImg1,1);
%     CmpImg1(1,:)=0;
%     CmpImg2=conj(flipud(CmpImg1(2:end,:)));
%     DblImg=[CmpImg1; CmpImg2];
%     DblSpm=ifft(DblImg);
%     fid = fopen(['Frame', FileName(end-6:end-4),'.dat'],'w');
%     fwrite(fid, DblSpm, 'float');fclose(fid);
    save(['Frame', FileName(end-6:end-4),'.mat'],'CaliedData');
end


%% save the frame diff image for further angiographic study
% IterationNum=4;
% SurfaceDepth=1;
% PenatratedDepth=600;
% Thickness=PenatratedDepth-SurfaceDepth+1;
% EnhancedIntensityImage=IntensityImage(SurfaceDepth:PenatratedDepth,:);%this one has problem because only one time cropping is needed.
% ScanNum=size(EnhancedIntensityImage,2)/IterationNum;
% % if get(handles.AngioByStd,'Value') && (FileName(end-4)=='B')
% StackImage=zeros(Thickness,ScanNum-ScanNum/10,IterationNum);
% TempImage=abs(EnhancedIntensityImage(:,ScanNum/10+1:ScanNum));%pls input the
% StackImage(:,:,1)=TempImage;
% for i=1:IterationNum-1
%     SubImage=EnhancedIntensityImage(:,ScanNum*(i+0.1)+1:ScanNum*(i+1));
%     ImageCC=normxcorr4smallshift(TempImage,SubImage,5, 5);
%     [MaxCC,Index]=max(ImageCC(:));
%     [RowShift,ColShift]=ind2sub(size(ImageCC),Index(1));
%     RowShift=RowShift-6;
%     ColShift=ColShift-6;
%     StackImage(:,:,i+1)=circshift(SubImage,[RowShift,ColShift]);
%     StackImage(:,:,i+1)=SubImage;
% end
% EnhancedIntensityImage=std(StackImage,1,3)./mean(StackImage,3);
%     
%% Filtering and resize and Ranging
if Procedure<=Const('Filter');
    EnhancedIntensityImage=IntensityImage(uint16(ImgPar.ImageFOV(1)*end+1:ImgPar.ImageFOV(2)*end),uint16(ImgPar.ImageFOV(3)*end+1:ImgPar.ImageFOV(4)*end));
%     EnhancedIntensityImage=LateralFreqFilter(EnhancedIntensityImage,IntPar.CutoffFreq,IntPar.CutoffAmp); 
    EnhancedIntensityImage=FiltImage(EnhancedIntensityImage, ImgPar.ImgFilterSize, ImgPar.ImgFilterType);    
    EnhancedIntensityImage=imresize(EnhancedIntensityImage,ImgPar.ImgResize.*size(EnhancedIntensityImage),'box');
end

% %% get PSF for different depth
% if debug
%     if isempty(AlineProfile)
%         AlineProfile=zeros(size(EnhancedIntensityImage,1),21);
%     end
% 
%     i=str2num(FileName(end-5:end-4));
%     AlineProfile(:,i)=mean(EnhancedIntensityImage,2);
% end    


%% display and save, colormap can be add here
if Procedure<=Const('Display');
    DisplayRange=GetDisplayRange(EnhancedIntensityImage, IntPar.IntDispRange);
    if ishandle(hIntensityImage)
% %         set(hIntensityImage,'CData', EnhancedIntensityImage); 
% %         set(hIntensityImage,'DisplayRange',DisplayRange);
%         set(0,'CurrentFigure',hIntensityImage);
%         imshow(EnhancedIntensityImage,'DisplayRange',DisplayRange);
%     else
%         hIntensityImage=figure('Name','IntensityImage','NumberTitle','off');
%         imshow(EnhancedIntensityImage,'DisplayRange',DisplayRange);impixelinfo;colormap gray;
%     end   
%         set(0,hIntensityImage);
        set(hIntensityImage,'CData',EnhancedIntensityImage); 
    else
        figure('Name','IntensityImage','NumberTitle','off');
        hIntensityImage=imshow(EnhancedIntensityImage,'DisplayRange',DisplayRange);
        impixelinfo;
    end

%    figure(2);hold on;plot(mean(EnhancedIntensityImage(:, 250:260)'));
    
    %% save
    if IntPar.AutoSaveFlt
        fid = fopen([FilePath,'OCTImage',FileName(end-6:end-4),'.dat'],'w');
        fwrite(fid, EnhancedIntensityImage, 'float');fclose(fid);
    end
    if IntPar.AutoSaveImg
        imwrite(uint16((EnhancedIntensityImage-DisplayRange(1))/(DisplayRange(2)-DisplayRange(1))*65535),[FilePath,'OCTImage',FileName(end-6:end-4),'.tif']);
    end
end