%caculate the dispersion as degree/mm for each wavelength
% close all;

depth = 7400
fftnum = 8192
pixelnum = 1376
n = 1;
thickness = 25;
iDepth = (0:depth/(fftnum-1):depth)-depth/2

sp = 1 .*gausswin(pixelnum, 3)';
% sp = sp0
c = 3*10^14;
wc = c/1.0*2*pi;
w0 = c/0.85*2*pi;
w1 = c/1.15*2*pi;
% sp = gausswin(pixelnum, 3)';
% sp = fliplr(sp)
d_w = (w0-w1)/(pixelnum-1);
w = c/1.15*2*pi:d_w:c/0.85*2*pi;
lambda = fliplr(c*2 *pi./w);
k1 = -0.0013267;% not sure about this value
k2 = -0.009;
k3 = -0.11;


    
% figure(1); plot(lambda, zeros(1,pixelnum));
% figure(2); plot(lambda, real(sp'));
% figure(3); plot((0:depth/(fftnum-1):depth)-depth/2 ,fftshift(abs(fft(sp', fftnum))))

d_phase = zeros(n+1, pixelnum);
spd = zeros(n+1, pixelnum);
for i = 1:n
    l = i*thickness*2*10^3;
    d_phase(i+1,:) =k1* 10^(-15)*(w-wc)*l + k2* 10^(-30)*(w-wc).^2 * l /2 + k3* 10^(-45)*(w-wc).^3 * l /6;    
    spd(i+1,:) = fliplr(sp).*exp(1j*(d_phase(i+1,:)));
end
d_phase = fliplr(d_phase);
figure;plot(x, d_phase');
spd = fliplr(spd);
psf = abs(fft(spd',fftnum)+1);
psf = fftshift(psf,1);
psf = psf /max(psf(:));
% psf = 20*log10(psf);
for i = 0 : n
    [m, im] = max(psf(:,i+1));
    i1 = interp1(psf(1:im,i+1), iDepth(1:im), m/2);
    i2 = interp1( psf(im:end,i+1),iDepth(im:end), m/2);
    FWHM(i+1) = i2-i1+1;
end
figure(1);
plot(FWHM/min(FWHM))
% figure(1);hold on; plot(lambda, d_phase);legend('0mm','10mm','20mm')%,'15mm','20mm','25mm');
% ylabel('rad'); xlabel('wavelength, um')
% figure(2);hold on; plot(lambda, real(spd'));legend('0mm','10mm','20mm')%','15mm','20mm','25mm');
% ylabel('intensity'); xlabel('wavelength, um')
% figure(3);hold on; plot((0:depth/(fftnum-1):depth)-depth/2, fftshift(psf,1));legend('0mm','10mm','20mm');
% ylabel('intensity'); xlabel('depth, um')
% 
figure(4); plot( (0:depth/(fftnum-1):depth)-depth/2, psf);
ylabel('intensity'); xlabel('depth, um')
legend('without window function','with window function')


ylabel('PSF broading ratio'); xlabel('Water thickness, mm')
