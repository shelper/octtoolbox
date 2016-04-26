function DeConjugatedData=DeConjugate(handles, CaliedData,DispersionPhaseError)

% DenoisedData=CaliedData;
%% need to redefine the filter
AlineNum=size(CaliedData,2);

DeConjWindow=get(handles.DeConjWindow,'String');
DeConjWindow=str2num(DeConjWindow(11:end));
ZeroPad=round(AlineNum*(1-diff(DeConjWindow))/diff(DeConjWindow));
Filter1=gausswin(AlineNum/2,3.6)';Filter1(end/2:end)=1;Filter2=fliplr(Filter1);
Filter=[Filter2,zeros(1,ZeroPad),Filter1];
Filter=imresize(Filter, [1,AlineNum]);
Filter=circshift(Filter,[1,round(mean(DeConjWindow)*AlineNum)]);
% Filter=zeros(1,AlineNum);
% Filter(1:end/2)=1;
% figure;plot(Filter);

fftData=ifft(CaliedData,[],2)*sparse(diag(Filter));
PosDataWithNeg=fft(fftData,[],2);
AccPosData=sparse(diag(exp(1i*DispersionPhaseError)))*PosDataWithNeg;
fftData=ifft(CaliedData,[],2)*sparse(diag(1-Filter));
NegDataWithPos=fft(fftData,[],2);
AccNegData=sparse(diag(exp(-1i*DispersionPhaseError)))*NegDataWithPos;

PosDataWithNeg=sparse(diag(exp(-1i*DispersionPhaseError)))*PosDataWithNeg;
NegDataWithPos=sparse(diag(exp(1i*DispersionPhaseError)))*NegDataWithPos;
clear CaliedData fftData;

DeConjItNum=str2double(get(handles.DeConjItNum,'String'));
PeakWidth=str2double(get(handles.PeakWidth,'String'));
StdThresh=str2double(get(handles.StdThresh,'String'));
PeakThresh=str2double(get(handles.PeakThresh,'String'));

for n=1:DeConjItNum
    %% added
    if n==1
        PosImgWithNeg=fftshift(ifft(PosDataWithNeg),1);
        SortedAbsImg=sort(abs(PosImgWithNeg),1); BGThreshold=median(SortedAbsImg(end/2,:));clear SortedAbsImg;        
        PosImg=FindOCTPeaks(PosImgWithNeg, PeakWidth, StdThresh, BGThreshold+0.1, PeakThresh);
        NegImgWithPos=fftshift(ifft(NegDataWithPos),1);
        NegImg=FindOCTPeaks(NegImgWithPos, PeakWidth, StdThresh, BGThreshold+0.1, PeakThresh);
%         AccPosImg=AccPosImg-NegImg;
%         AccNegImg=AccNegImg-PosImg;
    else
        NegDataWithPos=fft(ifftshift(PosImgWithNeg-PosImg,1));   
        PosDataWithNeg=fft(ifftshift(NegImgWithPos-NegImg,1));  
        
        NegDataWithPos=sparse(diag(exp(1i*DispersionPhaseError*2)))*NegDataWithPos;
        NegImgWithPos=fftshift(ifft(NegDataWithPos),1);
        NegImg=FindOCTPeaks(NegImgWithPos, PeakWidth, StdThresh, BGThreshold-0.1, PeakThresh+0.08);


        PosDataWithNeg=sparse(diag(exp(-1i*DispersionPhaseError*2)))*PosDataWithNeg;
        PosImgWithNeg=fftshift(ifft(PosDataWithNeg),1);        
        PosImg=FindOCTPeaks(PosImgWithNeg, PeakWidth, StdThresh, BGThreshold-0.1, PeakThresh+0.08);
%         figure;imshow(abs(PosImg)/100);
        
        AccPosData=AccPosData-fft(ifftshift(NegImg,1));        
        AccNegData=AccNegData-fft(ifftshift(PosImg,1));
    end    
    
%     PosDataWithNeg=fft(ifftshift(NegImgWithPos-NegImg,1));  
%     PosDataWithNeg=sparse(diag(exp(1i*DispersionPhaseError*2)))*PosDataWithNeg;
end
clear PosImg NegImg PosImgWithNeg NegImgWithPos PosDataWithNeg NegDataWithPos;
% PosImgWithNeg=fftshift(ifft(PosDataWithNeg),1);
% PosImg=FindOCTPeaks(PosImgWithNeg, PeakWidth, StdThresh, BGThreshold, PeakThresh);
% AccPosImg=AccPosImg+PosImg;   
DeConjugatedData=sparse(diag(exp(1i*DispersionPhaseError*2)))*AccNegData+conj(sparse(diag(exp(-1i*DispersionPhaseError*2)))*AccPosData);
DeConjugatedData=DeConjugatedData/2;
% AccNegData=fft(ifftshift(AccNegImg,1));
% AccPosData=fft(ifftshift(AccPosImg,1));
% DeConjugatedData=AccPosData+conj(AccNegData);
% DeConjugatedData=CaliedData-AccPosData;

% AccNegData=sparse(diag(exp(1i*DispersionPhaseError*2)))*AccNegData;
% % CaliedData=sparse(diag(exp(-1i*DispersionPhaseError)))*CaliedData;

% AccPosData=fft(ifftshift(AccPosImg,1));
% AccPosData=sparse(diag(exp(-1i*DispersionPhaseError*2)))*AccPosData;
% CaliedData=sparse(diag(exp(-1i*DispersionPhaseError)))*CaliedData;
% DeConjugatedData=CaliedData-AccPosData;
end

function OutputImage=FindOCTPeaks(InputImage, PeakWidth, StdThreshold,BGThreshold, PeakThresh)
    AbsImg=abs(InputImage);
    SortedAbsImg=sort(AbsImg,1);
    Filter=fspecial('average',[PeakWidth,1]);
    
    PeakThresh=median(SortedAbsImg(round(PeakThresh*end),:));  %this is set to 0.9*end:end
    StdImg=sqrt((imfilter(AbsImg.^2, Filter)-(imfilter(AbsImg,Filter)).^2))./imfilter(AbsImg, Filter);
    Mask=(AbsImg>BGThreshold).*(AbsImg>PeakThresh).*(StdImg>StdThreshold);
    
%     StdThreshold=median(SortedAbsImg(round(PeakThresh*end:end),:))*StdThreshold;  %this is set to 0.9*end:end
%     StdImg=sqrt((imfilter(AbsImg.^2, Filter)-(imfilter(AbsImg,Filter)).^2));
%     Mask=(AbsImg>BGThreshold).*(StdImg>repmat(StdThreshold, size(StdImg,1),1));

    Mask=(imfilter(Mask,Filter)>0); %Mask=min(Mask,1);  %figure;imshow(Mask);
    OutputImage=InputImage.*Mask;
end


function DeConjugate2
% for n=1:IterationNum
%     PosImg=fftshift(ifft(PosData),1);
% %     figure;imshow(20*log10(abs(PosImg)),[70,110]);
%     
%     RestNegImg=PosImg-FindOCTPeaks(PosImg, 5, 0.45,[0.9,1]);
%     % can be optimized by histogram
%     % increase Stdthreshold, remove more noise. aggressive.
%     % increase PeakThresh(1), remove more noise.
%     % should  be more and more aggresive since conjugate mostly removed
%     RestNegData=fft(ifftshift(RestNegImg,1));    
%     RestNegData=sparse(diag(exp(-1i*DispersionPhaseError*2)))*RestNegData;
%     RestNegImg=fftshift(ifft(RestNegData),1);   
% %     figure;imshow(abs(RestNegImg)/200000);impixelinfo;
% %     figure;imshow(20*log10(abs(RestNegImg)),[70,110]);
%     RestNegImg=FindOCTPeaks(RestNegImg, 5, 0.45,[0.9,1]);
%     % can be optimized by histogram   
%     % Decrease Stdthreshold, remove more noise. aggressive.
%     % Decrease PeakThresh(1), remove more noise.
%     % should  be more and more concervertive to avoid removal of true
%     % signal, thus, should increase these two numbers
%     RestNegData=fft(ifftshift(RestNegImg,1));  
%     RestNegData=sparse(diag(exp(1i*DispersionPhaseError*2)))*RestNegData;
%     
% %     AccuRestNegData=AccuRestNegData+RestNegData;
% %     PosImg=PosImg-ifft(RestNegData);
%     PosData=PosData-RestNegData;
% %     figure;imshow(abs(PosImg)/200000);impixelinfo;
% %     figure;imshow(fftshift(abs(PosImg)/100000,1));colormap hot;impixelinfo;    
% %     figure;imshow(abs(DenoisedNegImg)/200000);
% %     DenoisedPosData=fft(DenoisedPosImg);
% end
%     DeConjugatedData=PosData;
% end
end

function DeConjugate1
% DeConjugatedData=CaliedData-PosData;
% fftData=ifft(CaliedData,[],2);
% fftData(:,2:DeconjThresh*end)=0;
% DenoisedPosData=fft(fftData,[],2);% the filtering in lateral direction only work for one time.iteration wont help this process.
% DenoisedPosData=sparse(diag(exp(1i*DispersionPhaseError)))*DenoisedPosData;
% DenoisedPosImg=ifft(DenoisedPosData);
% DenoisedPosImg=DenoisedPosImg.*(abs(DenoisedPosImg)>1000);
% % DenoisedPosData=fft(DenoisedPosImg);
% 
% CaliedData=sparse(diag(exp(1i*DispersionPhaseError)))*CaliedData;
% CaliedImg=ifft(CaliedData);
% 
% for i=1:IterationNum
%     DenoisedNegImg=CaliedImg-DenoisedPosImg;
%     DenoisedNegData=fft(DenoisedNegImg);    
%     DenoisedNegData=sparse(diag(exp(-1i*DispersionPhaseError*2)))*DenoisedNegData;
%     DenoisedNegImg=ifft(DenoisedNegData);    
% %     figure;imshow(abs(DenoisedNegImg)/200000);
%     
%     DenoisedPosImg = circshift(ifftshift(flipud(conj(fftshift(DenoisedNegImg,1))),1),1); 
%     DenoisedPosImg=DenoisedPosImg.*(abs(DenoisedPosImg)>(500/i));
%     figure;imshow(abs(DenoisedPosImg)/200000);pause;
% %     DenoisedPosData=fft(DenoisedPosImg);
% 
% end
% 
% DenoisedPosImg=ifft(DenoisedPosData);
% DenoisedPosData=fft(DenoisedPosImg);

%     PosDispCorrNegImg=PosDispCorrImg-PosDispCorrPosImg;
%     PosDispCorrNegSig=fft(PosDispCorrNegImg);
%     NegDispCorrNegSig=PosDispCorrNegSig.*conj(DispSig.^2);    
%     NegDispCorrNegImg=ifft(NegDispCorrNegSig);    %step 2
%     
%     PosDispCorrPosImg=ifftshift(fliplr(conj(fftshift(NegDispCorrNegImg))));   %step 3
%     PosDispCorrPosImg=PosDispCorrPosImg.*(abs(PosDispCorrPosImg)>(200/i));
% %     PosDispCorrNegImg= PosDispCorrImg - PosDispCorrPosImg;  %step 4
%     plot(fftshift(abs(NegDispCorrNegImg)));pause;
    
    
%     fftData=ifft(DenoisedData,[],2);
%     fftData(:,2:DeconjThresh*end)=0;
%     DenoisedData=fft(fftData,[],2); 
    
%     if mod(i,2)
%         fftData=ifft(DenoisedData,[],2);
%         fftData(:,2:DeconjThresh*end)=0;
%         DenoisedData=fft(fftData,[],2);
%         if i==IterationNum
%             DenoisedData=sparse(diag(exp(1i*DispersionPhaseError)))*DenoisedData;
%             return;
%         else
% %           -sparse(diag(exp(-1i*DispersionPhaseError)))*
%             DenoisedData=CaliedData-DenoisedData;     
%         end
%     else
% %           DenoisedData=sparse(diag(exp(-1i*DispersionPhaseError)))*DenoisedData;
%         fftData=ifft(DenoisedData,[],2);
%         fftData(:,(1-DeconjThresh)*end+2:end)=0;
%         DenoisedData=fft(fftData,[],2);
%         if i==IterationNum
%             DenoisedData=sparse(diag(exp(-1i*DispersionPhaseError)))*DenoisedData;
%             DenoisedData=conj(DenoisedData);
%             return
%         else
%             DenoisedData=CaliedData-DenoisedData;
%         end      
%     end    
% end
end

function FindOCTPeaks1
% %% find peak by max value
% function OutputImage=FindPeaks(InputImage, PeakNum, HalfPeakWidth)
%     OutputImage=zeros(size(InputImage));
%     for m=1:PeakNum
%         [~,tmpIndex]=max(abs(InputImage));
%         tmpIndex=tmpIndex+(0:size(InputImage,1):numel(InputImage)-1);
%         for n=-HalfPeakWidth:HalfPeakWidth
%             Index=max(tmpIndex,1-n);
%             Index=min(Index,numel(InputImage)-n);            
%             OutputImage(Index+n)=max(InputImage(Index+n),OutputImage(Index+n));
%             InputImage(Index+n)=0;
%         end
%     end
% end
end
