% Calculate the correlation to see if that matches theoretical prediction
FileNames = dir('*.raw');
Style = ['r', 'g', 'b', 'k'];
DepthSpan = 40;
DepthRange = 151:450;

for iFile = 1:size(FileNames, 1)
    OCTData = GetOCTData(FileNames(iFile).name, 'OCT2000');
%     OCTData = GetFreqComp(OCTData, [200, 600], 2);
    CplxImg = fft(OCTData);
    CplxImg = CplxImg(DepthRange, :);
    OCTImg = abs(CplxImg);
%     figure;imshow(OCTImg,[30,100]);
    %% correlation by Doppler Phase
    WeightImg = OCTImg(:, 1:end-1) + OCTImg(:, 2:end);
    WeightImg = WeightImg .* (WeightImg>500);
    DopImg = angle(conj(CplxImg(:, 1:end-1)) .* CplxImg(:, 2:end));
    DopMean = sum(DopImg .* WeightImg)./sum(WeightImg);
    TempCorrP = sqrt(sum(bsxfun(@minus, DopImg, DopMean).^2 .*WeightImg)./sum(WeightImg));
    TempCorrP(isnan(TempCorrP)) = [];
    figure(1);hold on;plot(TempCorrP,Style(iFile));title('TempCorrP');
    DopImgCorrP(iFile) = mean(TempCorrP);
%     DiffImg = exp(1i*DopImg);
%     RealImg = real(DiffImg);
%     ImagImg = real(DiffImg);
%     figure;imshow(DopImg,[]);colormap jet;
%     DiffMedian = median(real(DiffImg))+1i*median(imag(DiffImg));
%     TempCorrC = sqrt(sum(abs(bsxfun(@minus, DiffImg,DiffMedian)).^2 .*WeightImg)./sum(WeightImg));
%     TempCorrC(isnan(TempCorrC)) = [];
%     figure(2);hold on;plot(TempCorrC,Style(iFile));title('TempCorrC');
%     DopImgCorrC(iFile) = mean(TempCorrC);
%     figure;plot(TempCorr);
    %% correlation by raw data
%     for iLine = 1:size(OCTData,2)-1
%         [TempCorr(iLine), Lag(iLine)] = max(xcorr(OCTData(:,iLine), OCTData(:,iLine+1), 3,'coeff'));
%     end
%     figure(3);hold on;plot(TempCorr,Style(iFile));
%     RawCorr(iFile) = mean(TempCorr);
    %% correlation by logscal image
%     for iLine = 1:size(LogImg,2)-1
%         [TempCorr(iLine), Lag(iLine)] = max(xcorr(LogImg(:,iLine), LogImg(:,iLine+1), 3,'coeff'));
%     end
%     LogImgCorr(iFile) = mean(TempCorr); 
    %% correlation by linear scale complex image
%     for iLine = 1:size(CplxImg,2)-1
%         [TempCorr(iLine), Lag(iLine)] = max(xcorr(CplxImg(:,iLine), CplxImg(:,iLine+1), 3,'coeff'));
%     end
%     figure(3);hold on;plot(TempCorr,Style(iFile));
%     CplxImgCorr(iFile) = mean(TempCorr);
    %% correlation by linear scale image
%     for iLine = 1:size(OCTImg,2)-1
%         TempCorr(iLine) = GetCorr(OCTImg(:, iLine), OCTImg(:, iLine+1), DepthSpan);
%     end
%     TempCorr(isnan(TempCorr)) = [];
%     Mask = (TempCorr>(mean(TempCorr)-0.1)) & (TempCorr<(mean(TempCorr)+0.1));
%     OCTImgCorr(iFile) = mean(TempCorr(Mask)); 
end

figure;plot(DopImgCorrP);
StepSize = [5.9,10.4,17.6,23.4];
PsfWidth = 16;
Xcorr0 = exp(-(StepSize.^2/PsfWidth.^2));
hold on, plot(StepSize/10, 'o-');


