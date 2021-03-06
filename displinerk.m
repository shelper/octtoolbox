%% system specs
k0 = (0: 1/1023:1);
pad0 = 2048;
k0_2048 = (0: 1/(pad0-1):1);
ks = k0.^2 * 0.156 + k0 *0.844 ; %from oct 2000
ks_2048 = k0_2048.^2 * 0.156 + k0_2048 *0.844 ;

for zd = 1:100:501  % dispersion match at different depth  
    disp = (2*pi* zd*k0 * 0.156 - 2*pi* zd*k0.^2*0.156);
    disp_2048 = (2*pi* zd *k0_2048 * 0.156 - 2*pi* zd *k0_2048.^2*0.156);
    color = ['r','g','b','y','k'];
    for n = 1:5 % different depth     
        z = n*30+150;
        % zero-padding of interferogram in Fourier domain
        sig = cos(2*pi* z *ks + disp ) .* gausswin(1024)';
        sigfft = fft(sig);
        sigfft_2048 = [sigfft, zeros(1, pad0-1024)];
        sig_2048 = real(ifft(sigfft_2048));
        % interpolation and reconstruction
        sig_linear = spline(ks_2048, sig_2048, k0_2048);
        disp_linear = spline(ks_2048, disp_2048, k0_2048);
        sig_linear = sig_linear.*exp(-1i*disp_linear);
%         figure(zd); hold on;plot((abs(fft(sig_linear))));
        figure(1); hold on;plot(20*log10(abs(fft(sig_linear,2048))),'color',color(n));
    end
    xlim([0, 512]);
%     ylim([0, 500]);
    ylim([-60, 50]);
end