% GetDiffPhaseStd
function DopCorr = GetDiffPhaseStd(CplxImg, PhaseThreshMin, PhaseThreshMax) 
    OCTImg = abs(CplxImg);
%     figure;imshow(OCTImg,[30,100]);
    %% correlation by Doppler Phase
    WeightImg = OCTImg(:, 1:end-1) + OCTImg(:, 2:end);
    WeightImg = WeightImg .* (WeightImg>PhaseThreshMin);
    WeightSum = sum(WeightImg) + (sum(WeightImg) ==0);
    DopImg = angle(conj(CplxImg(:, 1:end-1)) .* CplxImg(:, 2:end));
    DopMean = sum(DopImg .* WeightImg)./WeightSum;
    DopCorr = sqrt(sum(bsxfun(@minus, DopImg, DopMean).^2 .*WeightImg)./WeightSum);
%     figure(1);hold on;plot(DopCorr,Style(iFile));title('TempCorrP');
%     DopImgCorrP(iFile) = mean(DopCorr);