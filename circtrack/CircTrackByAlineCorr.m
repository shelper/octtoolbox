%% get system parameters
delta_x = 0.010; %mm, x100um
aline_num = 1024;
frame_num = 50;
aline_rate = 27000; %lines/s
width_psf = 0.008; %mm, x1000um
vc = delta_x*aline_rate; %mm/s
DepthRange = 151:550;
BgCorr = zeros(1, aline_num-1);
MvCorr = zeros(1, aline_num-1);

%% Get the first frame as motion free background
FileName = ['.\0518\', 'OCTImg', sprintf('%.4d', 0), '.raw'];
OCTData = GetOCTData(FileName, 'OCT2000');
% OCTData = bsxfun(@times, OCTData,gausswin(2048,1));
CplxImg = fft(OCTData);
CplxImg = CplxImg(DepthRange, :);
%% Get Aline correlation for tracking
% CircData = abs(CplxImg);
% for iLine = 1:size(CircData,2)-1
%     BgCorr(iLine) = GetCorr(CircData(:, iLine), CircData(:, iLine+1), 40);
% end
% BgCorr(isnan(BgCorr)) = [];

%% Get Diff Phase Variation for tracking
BgCorr = GetDiffPhaseStd(CplxImg, 500);

am_calc = zeros(3, frame_num);
vm_calc = zeros(3, frame_num);
%% get following frames to calculate the motion
for i = 1:frame_num-1
    FileName = ['.\0518\', 'OCTImg', sprintf('%.4d', i), '.raw'];
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
    displace = MvCorr-BgCorr;
%     A = displace(1:aline_num/4);
%     A = mean(A(A>mean(A)-2*std(A) & A<mean(A)+2*std(A)));
%     B = displace(aline_num/4+1:aline_num/2);
%     B = mean(B(B>mean(B)-2*std(B) & B<mean(B)+2*std(B)));
%     C = displace(aline_num/2+1:aline_num*3/4);
%     C = mean(C(C>mean(C)-2*std(C) & C<mean(C)+2*std(C)));
%     D = displace(aline_num*3/4+1:end);
%     D = mean(D(D>mean(D)-2*std(D) & D<mean(D)+2*std(D)));
    A = mean(displace(1:aline_num/4));
    B = mean(displace(aline_num/4+1:aline_num/2));
    C = mean(displace(aline_num/2+1:aline_num*3/4));
    D = mean(displace(aline_num*3/4+1:end));
    %% use A, B
    am_calc(1, i) = atan(B/A);
    if A<0 && B<0
        am_calc(1, i) = am_calc(1, i)-pi ;
    elseif  A<0 && B>0 
        am_calc(1, i) = am_calc(1, i)+pi ;
    end
    vm_calc(1, i) = abs(A/(2*vc*sin(am_calc(1, i))));
%     am_calc(1, i) = angle(exp(1i*(am_calc(1, i)+pi/4)));
    %% use B, C
    am_calc(2, i) = atan(-B/C);
    if B<0 && C>0
        am_calc(2, i) = am_calc(2, i)-pi ;
    elseif  C>0 && B>0 
        am_calc(2, i) = am_calc(2, i)+pi ;
    end
    vm_calc(2, i) = abs(C/(2*vc*sin(am_calc(2, i))));
%     am_calc(2, i) = angle(exp(1i*(am_calc(2, i)+pi/4)));
    %% use C, D
    am_calc(3, i) = atan(D/C);
    if D>0 && C>0
        am_calc(3, i) = am_calc(3, i)-pi ;
    elseif  C>0 && D<0 
        am_calc(3, i) = am_calc(3, i)+pi ;
    end
    vm_calc(3, i) = abs(C/(2*vc*sin(am_calc(3, i))));
    figure(10); hold on; plot([A B C D]);
%     am_calc(3, i) = angle(exp(1i*(am_calc(3, i)+pi/4)));
%     BgCorr = MvCorr;
end  

figure;plot(vm_calc');
figure;plot(am_calc');
% figure;plot(exp(1i*am_calc(:)));
vm_calc(vm_calc>50) = 0;
figure;plot(cumsum(vm_calc.*exp(1i*(am_calc)),2)');
hold on;plot(cumsum(median(vm_calc).*exp(1i*median(am_calc))),'k');

%% calculate averaged displace of each Quadrant
%% calculate the speed and the displace


