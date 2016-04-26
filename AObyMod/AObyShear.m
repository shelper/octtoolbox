close all;
clear all;
arr_size  = 100;
phase_sam = randn(arr_size)*2*pi;
phase_sam = imfilter(phase_sam, fspecial('gaussian', 6,2) );
beam_profile = fspecial('gaussian', 100,30);
beam_profile =beam_profile / max(beam_profile(:));
sam = exp(1i*phase_sam).*beam_profile;
figure;imshow(angle(sam),[-pi, pi]);colormap jet;

shift = [1,0];
phase_shear = zeros(arr_size);
for m = 1:arr_size
    for n = 1:arr_size
        phase_shear(m, n) = (m+n) * pi/2;    
    end
end
sam_shear = circshift(sam, shift);
sam_shear = exp(1i*phase_shear).* sam_shear;

noise = randn(size(sam))*0.05;
xcc = (sam+sam_shear).*conj(sam+sam_shear) + noise;
xcc_fft = fft2(xcc);
xcc_filt = zeros(size(xcc_fft));
xcc_filt(51:end-10,51:end-10) = xcc_fft(51:end-10,51:end-10);
xcc_filt = ifft2(xcc_filt);
if shift(1) == 2
    phase_diff = angle(sqrt(xcc_filt.* exp(1i*phase_shear)));
else
    phase_diff = angle(xcc_filt.* exp(1i*phase_shear));
end
phase_diff_center = phase_diff.*(beam_profile>=0.2);

phase_diff_fft = fft(phase_diff);
phase_fft = [zeros(1,100);bsxfun(@rdivide, phase_diff_fft(2:end,:), 1i/(5*pi)*[1:50,-49:-1]')];
phase_recon = ifft(phase_fft);
figure;imshow(phase_recon,[-pi, pi]);colormap jet;

phase_error = angle(exp(1i*(phase_recon-angle(sam))));
phase_error = phase_error-mean(phase_error(:));
figure;imshow(phase_error,[-pi, pi]);colormap jet;

figure;imshow(fftshift(abs(fft2(sam))),[0,1000]);
sam_recon = exp(1i*phase_error).*beam_profile;
figure;imshow(fftshift(abs(fft2(sam_recon))),[0,1000]);

% phase_recon = angle(exp(1i*cumsum(phase_diff, 1)));
% figure;imshow(phase_recon,[]);colormap jet;


