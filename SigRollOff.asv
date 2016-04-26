N = 20480;
psf_width = 18;
psf_thresh = 1/exp(2);
alpha = N * sqrt(2*log(1/psf_thresh))/psf_width;
psf=gausswin(N, alpha);
psf_fall=20*log10(abs(fft(psf)));
psf_fall = psf_fall(1:end/20);
psf_fall = psf_fall - max(psf_fall);
figure;plot(psf_fall, 'r');

pix=zeros(20480,1);
pix(1:10) = 1;
pix_fall=20*log10(abs(fft(pix)));
pix_fall = pix_fall(1:end/20);
pix_fall = pix_fall - max(pix_fall);
hold on;plot(pix_fall,'g')

k0 = (0:1/2048:1);
ks = k0.^2 * 0.156 + k0*0.844 ; %from oct 2000
for depth = 1:1024
    sig1 = cos(2*pi* depth *ks);
    sig2 = cos(2*pi* depth *k0); 
    sig3 = spline(ks, sig1, k0);
    psf2 = 20*log10(max(abs(fft(sig2, 4096))));
    psf3 = 20*log10(max(abs(fft(sig3, 4096))));
    resamp_fall(depth) = psf3-psf2;
end
hold on;plot(resamp_fall,'k')

hold on;plot(pix_fall + psf_fall+resamp_fall','b')
xlim([0, 1024])

legend('psf','pixelization','resample', 'combined');