%% get system parameters
delta_x = 0.010; %mm, x100um
aline_num = 1024;
frame_num = 50;
aline_rate = 27000; %lines/s
width_psf = 0.008; %mm, x1000um
vc = delta_x*aline_rate; %mm/s
DepthRange = 151:550;

BgCorr = zeros(1, aline_num-1);
%% Get the first frame as motion free background
FileName = ['.\bg\', 'OCTImg', sprintf('%.4d', 49), '.raw'];
OCTData = GetOCTData(FileName, 'OCT2000');
CplxImg = fft(OCTData);
CplxImg = CplxImg(DepthRange, :);
BgCorr = BgCorr + GetDiffPhaseStd(CplxImg, 500);  
% figure;plot(BgCorr);

%% get following frames to calculate the motion
MvCorr = zeros(1, aline_num-1);
for n = 1:frame_num
    n
    FileName = ['.\move\', 'OCTImg', sprintf('%.4d', n-1), '.raw'];
    OCTData = GetOCTData(FileName, 'OCT2000');
%     OCTData = bsxfun(@times, OCTData,gausswin(2048,1));
    CplxImg = fft(OCTData);
    CplxImg = CplxImg(DepthRange, :);
%% Get Aline correlation for tracking
%     MovedData = abs(CplxImg);
%     for iLine = 1:size(MovedData,2)-1
%         MvCorr(iLine) = GetCorr(MovedData(:, iLine), MovedData(:, iLine+1), 40);
%     end
%     MvCorr(isnan(MvCorr)) = 0;
%     displace = (aline_rate^2) * real((width_psf *sqrt(log(1./MvCorr))).^2-(width_psf * sqrt(log(1./BgCorr))).^2) ;
%% Get Diff Phase Variation for tracking
    MvCorr = GetDiffPhaseStd(CplxImg, 500);
    CorrDiff = MvCorr-BgCorr;
    CorrDiffMask = abs(CorrDiff)<0.3;
    CorrDiff(~CorrDiffMask) = 0;
    A(n) = sum(CorrDiff(1:aline_num/4))/sum(CorrDiffMask(1:aline_num/4));
    B(n) = sum(CorrDiff(aline_num/4+1:aline_num/2))/sum(CorrDiffMask(aline_num/4+1:aline_num/2));
    C(n) = sum(CorrDiff(aline_num/2+1:aline_num*3/4))/sum(CorrDiffMask(aline_num/2+1:aline_num*3/4));
    D(n) = sum(CorrDiff(aline_num*3/4+1:end))/sum(CorrDiffMask(aline_num*3/4+1:end));
end     
MoveVel = ((D-B) + 1i * (C-A)) ;
MoveVel = MoveVel* exp(1i*3*pi/4);
figure;plot(MoveVel,'-+');
figure;plot(cumsum(MoveVel), '-+');

%% calculate averaged displace of each Quadrant
%% calculate the speed and the displace


