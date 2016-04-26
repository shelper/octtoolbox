function CaliedData=PreProcess(PrePar,RawData, FilePath, FileName)
global DispersionPhaseError CaliCoeff Ref FRHRatio ;
debug=0;
%% split Spectrum OCT
% M = 10;
% SubWidth = floor(size(RawData,1)/M);
% RawData = RawData(end - SubWidth:end, :);

[PixelNum, AlineNum]=size(RawData);
% Ref = 0;
% if isempty(Ref)
%     Ref = mean(RawData,2);
% end
% RawData=RawData-repmat(Ref, 1, AlineNum);    

% RawData = RawData';
% spr = mean(RawData);
% sps = mean(abs(RawData-repmat(spr,AlineNum,1)));
% sps = sps / max(sps);
% RawData = RawData';
% RawData(1622,:) = (RawData(1621,:) + RawData(1623,:))/2;
% RawData=RawData.*repmat(barthannwin(PixelNum),1,AlineNum);

% if isempty(Ref)
%     Ref = mean(RawData, 2);
%     Ref = repmat(Ref, 1, size(RawData,2));
% end
% CaliedData = RawData - Ref; 
% CaliedData=RawData;
CaliedData=RemoveRef(RawData,PrePar.ReferenceMode);

% CaliCoeff=fminsearch(@(CaliCoeff) GetImageEntropy4Cali(CaliedData,CaliCoeff),[0,-0]);
% CaliCoeff = [CaliCoeff, 1-CaliCoeff(1),-CaliCoeff(2)]
% CaliK=(-1:2/2047:1)';
% OrigK=polyval(CaliCoeff,CaliK);
% CaliedData=interp1( OrigK,RawData, CaliK, 'linear', 'extrap');

%% calibration
% if ~isempty(CaliCoeff)
%     CaliedData=sparse(diag(CaliCoeff(2,:)))*RawData(CaliCoeff(1,:),:)...
%                +sparse(diag(CaliCoeff(4,:)))*RawData(CaliCoeff(3,:),:);
% else
%     CaliedData=RawData;
% end

%% load dispersion and handling
if PrePar.GenDispCoeff
%     DispCoeff=[-300,0];
    DispCoeff = fminsearch(@(DispCoeff) GetImageEntropy(CaliedData,DispCoeff),[-300, 0]); 
    DispersionPhaseError=GetNonLinearPhase(DispCoeff,PixelNum);
    save(['C:\Users\zyuan\Dropbox\research\OCTToolbox\Dispersion\',PrePar.DispFileName],'DispersionPhaseError');
elseif ~isempty(PrePar.DispFileName)
    load(['C:\Users\zyuan\Dropbox\research\OCTToolbox\Dispersion\',PrePar.DispFileName]);
else
    DispersionPhaseError=[];
end

% DispersionPhaseError = -GenDispCoeff(sum(CaliedData,2),[110, 170]);
% DispersionPhaseError=GetNonLinearPhase([136, 0],2048);

if ~isempty(DispersionPhaseError)
%     DispersionPhaseError=DispersionPhaseError(4:end-1);
    CaliedData = sparse(diag(exp(1i*DispersionPhaseError)))*CaliedData;
end

% CaliedData=RemoveRef(CaliedData,PrePar.ReferenceMode);
%% spectroscopic OCT     
%     RawData(1:end/2,:)=0;
%     RawData(end/2+1:end,:)=RawData(end/2+1:end,:).*repmat(barthannwin(size(RawData,1)/2),1,size(RawData,2));
%     RawData(end/2+1:end,:)=0;
%     RawData(1:end/2,:)=RawData(1:end/2,:).*repmat(barthannwin(size(RawData,1)/2),1,size(RawData,2));
% if PrePar.CorrectPixelShift
%     RawData=CorrectPixelShift(RawData,10);
% end
%% Deconjugate by LPM
% if PrePar.DeConjugate
%     CaliedData=DeConjugateByLPM(PrePar, CaliedData);
% end
%% phaase jittering correction
% if strcmpi(PrePar.SystemID(1:5), 'ssOCT') 
%         CaliedData=CorrectPixelShift(CaliedData,4,FRHRatio);
% end
% CaliedData=RemoveRef(CaliedData,PrePar.ReferenceMode);

if debug
    fid = fopen(['OCTData',FileName(end-6:end-4),'.dat'],'w');
    fwrite(fid, CaliedData, 'float');fclose(fid);
end
% CaliedData=CorrectPixelShift(CaliedData,4,FRHRatio);

%     fid = fopen([FilePath,'FRHRatio.dat']);
%     if fid ~=-1
%         FRHRatio = fread(fid,'double');
%         fclose(fid);
%     else
%         FRHRatio = ones(1376,1);
%     end

%% Deconjugate by Bes   
if strcmpi(PrePar.DeConjugate,'DeConjByBes') 
    CaliedData = CorrectPixelShift(CaliedData,2,1);
    CaliedData=DeConjugateByBes(CaliedData,PrePar.DeConjOrder,FilePath);
elseif strcmpi(PrePar.DeConjugate,'DeConjByLPM')
    CaliedData=DeConjugateByLPM(PrePar, CaliedData);
end

% CaliedData=CaliedData(353:end,:);
% CaliedData=CaliedData.*repmat(barthannwin(size(CaliedData,1)),1,size(CaliedData,2));

%% reshape the spectrum
% CaliedData=CaliedData.*repmat(barthannwin(size(CaliedData,1))+0.001,1,size(CaliedData,2));

%% debug scripts
% if debug
%     CaliedData=repmat(circshift(gausswin(size(CaliedData,1)),300),[1,size(CaliedData,2)]).*CaliedData;%Spectrum Reshaping
%     CaliedData=repmat(gausswin(size(CaliedData,1)),[1,size(CaliedData,2)]).*CaliedData;%Spectrum Reshaping
%     CaliedData=GetFreqComp(CaliedData,[190 220]);
%     foo=angle((fft(CaliedData(:,510:520),2048)));
% 
%     foo=foo./repmat(max(foo),[size(foo,1),1]);
%     figure;plot(foo);
% 
%     % foo=foo([1:100,638:738,1276:1376],:);
%     foo=unwrap(angle(hilbert(foo)));
%     % figure;plot(foo(:,500:503));
%     PhaseDiff=(sum(foo(:,1:2:end),2)-sum(foo(:,2:2:end),2))/512;
%     % PhaseDiff=PhaseDiff-PhaseDiff(1);
%     if sum(PhaseDiff)<0
%         PhaseDiff=-PhaseDiff;
%     end
%     % 
%     figure;plot(PhaseDiff);
% end

%% Deconjugate by Dispersion
% if PrePar.DEFR
%     CaliedData=DeConjugateByDisp(PrePar, CaliedData,DispersionPhaseError);
% end

end

function CaliedData=DeConjugateByLPM(PrePar, CaliedData)
    %% need to redefine the filter
    AlineNum=size(CaliedData,2);
    ZeroPad=round(AlineNum*(1-diff(PrePar.DeConjWindow))/diff(PrePar.DeConjWindow));
    Filter1=gausswin(AlineNum/2,3.6)';Filter1(end/2:end)=1;Filter2=fliplr(Filter1);
    Filter=[Filter2,zeros(1,ZeroPad),Filter1];
    Filter=imresize(Filter, [1,AlineNum]);
    Filter=circshift(Filter,[1,round(mean(PrePar.DeConjWindow)*AlineNum)]);
    fftData=fft(CaliedData,[],2);
    CaliedData=ifft(fftData*sparse(diag(Filter)),[],2);
end

function CaliedData=DeConjugateByDisp(PrePar, CaliedData,DispersionPhaseError)
    %% need to redefine the filter
    for n=1:Const('DEFRItNum')
        PosImgWithNeg=fftshift(ifft(CaliedData),1);
        SortedAbsImg=sort(abs(PosImgWithNeg),1); BGThreshold=SortedAbsImg(end/2,:);clear SortedAbsImg;        
        PosImg=FindOCTPeaks(PosImgWithNeg, PrePar.PeakWidth, PrePar.StdThresh, BGThreshold+0.1, PrePar.PeakThresh);

        NegDataWithPos=fft(ifftshift(PosImgWithNeg-PosImg,1));           
        NegDataWithPos=sparse(diag(exp(1i*DispersionPhaseError*2)))*NegDataWithPos;
        NegImgWithPos=fftshift(ifft(NegDataWithPos),1);
        NegImg=FindOCTPeaks(NegImgWithPos, PrePar.PeakWidth, PrePar.StdThresh, BGThreshold-0.1, PrePar.PeakThresh+0.08);
        NegData=fft(ifftshift(NegImg,1));
        NegData=sparse(diag(exp(1i*DispersionPhaseError*2)))*NegData;

        CaliedData=CaliedData-NegData;        
    end    
end

function Entropy = GetImageEntropy(DispedData, DispCoeff)
    DispersionPhaseError=GetNonLinearPhase(DispCoeff,size(DispedData,1));
    DispedData=sparse(diag(exp(1i*DispersionPhaseError)))*DispedData;
    Image = abs(ifft(DispedData));
    Image=Image(1:end/2-100,:);
    Image=Image/sum(Image(:));
    Entropy = Image.*log(Image);
    Entropy = -sum(Entropy(:));
end

function CaliedData=CorrectDisp(CaliedData, DispersionPhaseError)
    Image = abs(ifft(CaliedData));
    Entropy1=sum(Image(:));
    NonDispedData=sparse(diag(exp(1i*DispersionPhaseError)))*CaliedData;
    Image = abs(ifft(NonDispedData));
    Entropy2=sum(Image(:));
    if Entropy1>Entropy2
        CaliedData = NonDispedData;
    else
        CaliedData = sparse(diag(exp(-1i*DispersionPhaseError)))*CaliedData;
    end
end

function NonLinearPhase = GetNonLinearPhase(DispCoeff,SpectrumWidth)
    WaveNumber = ((1:SpectrumWidth)-SpectrumWidth/2)/SpectrumWidth;
    DispPhase=DispCoeff(1)*WaveNumber.^2+DispCoeff(2)*WaveNumber.^3;
    LinearPhase=DispPhase(1):(DispPhase(end)-DispPhase(1))/(SpectrumWidth-1):DispPhase(end);
    if isempty(LinearPhase)
        LinearPhase=0;
    end
    NonLinearPhase=DispPhase-LinearPhase;
end

function OutputImage=FindOCTPeaks(InputImage, PeakWidth, StdThreshold,BGThreshold, PeakThresh)
    AbsImg=abs(InputImage);
    SortedAbsImg=sort(AbsImg,1);
    Filter=fspecial('average',[PeakWidth,1]);
    
    StdImg=sqrt((imfilter(AbsImg.^2, Filter)-(imfilter(AbsImg,Filter)).^2))./imfilter(AbsImg, Filter);
%     PeakThresh=median(SortedAbsImg(round(PeakThresh*end),:));  %this is set to 0.9*end:end
%     Mask=(AbsImg>BGThreshold).*(AbsImg>PeakThresh).*(StdImg>StdThreshold);
    PeakThresh=repmat(SortedAbsImg(round(PeakThresh*end),:), size(SortedAbsImg,1),1);  %this is set to 0.9*end:end
    Mask=(AbsImg>PeakThresh).*(StdImg>StdThreshold);
    
%     StdThreshold=median(SortedAbsImg(round(PeakThresh*end:end),:))*StdThreshold;  %this is set to 0.9*end:end
%     StdImg=sqrt((imfilter(AbsImg.^2, Filter)-(imfilter(AbsImg,Filter)).^2));
%     Mask=(AbsImg>BGThreshold).*(StdImg>repmat(StdThreshold, size(StdImg,1),1));

    Mask=(imfilter(Mask,Filter)>0); %Mask=min(Mask,1);  %figure;imshow(Mask);
    OutputImage=InputImage.*Mask;
end