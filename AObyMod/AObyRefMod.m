close all;
arr_size  = 100;
phase_sam = randn(arr_size)*2*pi;
phase_sam = imfilter(phase_sam, fspecial('gaussian', 6,2) );
beam_profile = fspecial('gaussian', 100, 25);
beam_profile =beam_profile / max(beam_profile(:));
sam = exp(1i*phase_sam).* beam_profile;
figure;imshow(angle(sam),[]);colormap jet;

phase_ref = zeros(arr_size);
for m = 1:arr_size
    for n = 1:arr_size
        phase_ref(m, n) = (m+n) * pi/2;
    end
end
ref = exp(1i*phase_ref).* beam_profile;

xcc = (sam+ref).*conj(sam+ref);
% figure;imshow((xcc), []);impixelinfo
% xcc = imrotate(xcc,45,'crop');
% figure;imshow((xcc), []);impixelinfo

xcc_fft = fft2(xcc);
figure;imshow(abs(xcc_fft), []);impixelinfo
xcc_filt = zeros(size(xcc_fft));
xcc_filt(5:50,5:50) = xcc_fft(5:50,5:50);

% xcc_filt(51-20:51+20,51-20:51+20) = xcc_fft(51-20:51+20,51-20:51+20);
xcc_filt = ifft2(xcc_filt);
% phase_ref = imrotate(phase_ref,45,'crop');
xcc_filt = xcc_filt.* exp(-1i*phase_ref);
% xcc_filt = imrotate(xcc_filt,-45,'crop');

figure;imshow(-angle(xcc_filt),[]);colormap jet;
