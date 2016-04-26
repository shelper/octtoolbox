a= gausswin (700, 3)*1.5;
b= gausswin(1000);
c=b;
c(1:700) = c(1:700) +a(1:700);
c(1:700) = c(1:700) -a(700);
c=c/max(c);
d = sqrt(flipud(c).*c);
d = d/max(d)
figure;plot(c,'b', 'LineWidth', 4)
hold on;plot(flipud(c),'r', 'LineWidth', 4)
hold on;plot(d,'--k', 'LineWidth', 4)


[~,i1]= min(abs(wabs0(:, 1)- 900));
[~,i2]= min(abs(wabs0(:, 1)- 1100));
wabs = wabs0(i1:i2,2);
lambda = 900:0.1:1100;
k0=1./lambda;
k = k0(1):(k0(end)-k0(1))/2000:k0(end);
wabs = interp1(k0, wabs, k);
fftnum = 4096;
res = 1/(k(1)-k(end))/2000;
res = res*size(wabs,2)/fftnum;
dr = 40;
depth = ((0:dr+2) -dr/2-1)* res;
for i = 0.945:0.3:3
    lsp = gausswin(size(wabs,2), i)';
    h_sp = figure();
    figure(h_sp);subplot(2,1,1);plot(lambda, lsp/max(lsp),'k', 'LineWidth', 4);
    %psf of light source
    lpsf = fftshift(abs(fft(lsp,fftnum)));
    [v, j]=max(lpsf);
    j1 = interp1(lpsf(j-5:j), j-5:j, v/2);j2 = interp1(lpsf(j:j+5), j:j+5, v/2);
    fwhm1 = (j2-j1) * res;
    lpsf = 20*log10(lpsf(fftnum/2-dr/2:fftnum/2+dr/2+2)/max(lpsf));    
    figure(h_sp);subplot(2,1,2);plot(depth, lpsf,'k', 'LineWidth', 2);
    %cross spectrum wo reshaping
    ssp = wabs.*lsp;
    figure(h_sp); subplot(2,1,1);hold on;plot(lambda,ssp/max(ssp), 'r', 'LineWidth', 4);
    csp = sqrt(ssp .* lsp); 
    figure(h_sp); subplot(2,1,1);hold on;plot(lambda,csp/max(csp), 'b', 'LineWidth', 4); 
    %psf of cross spectrum wo reshaping
    cpsf = fftshift(abs(fft(csp,fftnum)));
    [v, j]=max(cpsf);
    j1 = interp1(cpsf(j-5:j), j-5:j, v/2);j2 = interp1(cpsf(j:j+5), j:j+5, v/2);
    fwhm2 = (j2-j1) * res;
    cpsf = 20*log10(cpsf(fftnum/2-dr/2:fftnum/2+dr/2+2)/max(cpsf));    
    figure(h_sp); subplot(2,1,2);hold on;plot(depth, cpsf, 'b', 'LineWidth', 2);    
    %reshaped cross spectrum
    csp = gausswin(size(wabs, 2), 2.145)';
    figure(h_sp); subplot(2,1,1);hold on;plot(lambda,csp/max(csp),'--k', 'LineWidth', 4);
    %psf of cross spectrum with reshaping
    cpsf = fftshift(abs(fft(csp,fftnum)));
    [v, j]=max(cpsf);
    j1 = interp1(cpsf(j-5:j), j-5:j, v/2);j2 = interp1(cpsf(j:j+5), j:j+5, v/2);
    fwhm3 = (j2-j1) * res;
    cpsf = 20*log10(cpsf(fftnum/2-dr/2:fftnum/2+dr/2+2)/max(cpsf));    
    figure(h_sp);subplot(2,1,2); hold on;plot(depth, cpsf, '--k', 'LineWidth', 2);  
    %set the range of the plots
    subplot(2,1,1);
    ylim([0, 1.05]);
    subplot(2,1,2);
    xlim([min(depth),max(depth)]);
    ylim([-80, 0]);
    %reshaped reference
    rsp = csp.^2./ssp;
    figure(h_sp); subplot(2,1,1);hold on;plot(lambda,rsp/max(rsp),'g', 'LineWidth', 4);
    %set figure titles:
    subplot(2,1,1); title('spectral shape');
    xlabel('wavelength(nm)');ylabel('power(a.u.)');
    legend('source spectrum', 'sample spectrum', 'unreshaped spectrum', 'reshaped spectrum', 'reference spectrum');
    subplot(2,1,2); title('psf');
    xlabel('depth(um)');ylabel('intensity');
    legend(['psf of source, FWHM:', num2str(fwhm1), 'um'],  ... 
        ['unreshaped psf, FWHM:', num2str(fwhm2), 'um'], ...
        ['reshaped psf, FWHM:', num2str(fwhm3), 'um']);
end


% function fwhm = get_fwhm(curve)
% [v, i]=max(curve)/2;
% i1 = interp1(curve(i-5:i), i-5:i, v/2);i2 = interp1(curve(i:i+5), i:i+5, v/2);
% fwhm = (i2-i1) * res;















