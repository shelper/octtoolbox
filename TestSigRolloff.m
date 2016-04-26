close all;
AlineNum = 512;
PixelNumber = 2048;

fid1 = fopen('T:\ZYtemp\2012-12-17 SIGNAL ROLL OFF\0222\octimg0068.RAW');
% rawdata = fread(fid1, [PixelNumber , AlineNum], 'uint16');
% ks0 = load('T:\ZYtemp\2012-11-28 2048 camera evl\newcam.mat');   
% ks0 = ks0.phase';
% ref = median(rawdata, 2);

foo =  fread(fid1, 3, 'int');
rawdata = fread(fid1, [PixelNumber , AlineNum], 'uint16');
ks0 = fread(fid1, [PixelNumber , 1], 'float32')';
ref = fread(fid1, [PixelNumber , 1], 'double');

fclose(fid1);
rawdata = bsxfun(@minus, ref, rawdata);

% fid2 = fopen('C:\Topcon\Projects\2011-09-06 DOCT4THQ\RAW10\Rescale.clb');
% ks0= fread(fid2, PixelNumber, 'float32')';
% fclose(fid2);

%% linear
k0 = linspace(0,1,PixelNumber);
ks = polyval(polyfit(ks0, k0, 3), k0);
% ks = ks - min(ks);
% ks = ks/max(ks);

data_calied1=interp1( ks, rawdata, k0,'linear','extrap');% for new cam
% data_calied1=interp1( k0, rawdata, ks);% for oct2000
img_cpx =fft(data_calied1);
img = 20*log10(abs(img_cpx(100:1024,:)));
figure;imshow(img,[max(img(:))-45, max(img(:))]);
figure(10);hold on;plot(mean(img(:, 210:260),2),'r');

%% spline
data_calied2=interp1( ks, rawdata, k0, 'spline','extrap');% for new cam
data_calied2(2047:end,:)=0;
% data_calied2=interp1( k0, rawdata, ks, 'spline');% for oct2000
img_cpx =fft(data_calied2);
img = 20*log10(abs(img_cpx(100:1024,:)));
figure;imshow(img,[max(img(:))-45, max(img(:))]);
figure(10);hold on;plot(mean(img(:, 210:260), 2),'g');

%% KB-window
rsn = 3072;
k0 = linspace(0,1,rsn);
[win, ind, psf_corr] = CalibrateK(ks, k0, 5, 'besselwin'); % for oct2000

r1=rawdata(ind, :);
r2=rawdata(ind+1, :);
r3=rawdata(ind+2, :);
data_calied3=bsxfun(@times, r1,win(1,:)') +bsxfun(@times, r2,win(2,:)')+bsxfun(@times, r3,win(3,:)');

img_cpx =fft(data_calied3);
img = 20*log10(abs(img_cpx(100:1024,:)));
img = bsxfun(@plus, 20*log10(psf_corr(100:1024))', img);

figure;imshow(img,[max(img(:))-4, max(img(:))]);
figure(10);hold on;plot(mean(img(:, 210:260), 2),'b');

figure(10);legend('linear','spline', 'kb window');
