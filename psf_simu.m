lambda0 = 920;
lambda1 = 1180;
FWHM = 200;
[~,i1]= min(abs(wabs0(:, 1)- lambda0));
[~,i2]= min(abs(wabs0(:, 1)- lambda1));
wabs = wabs0(i1:i2,2)';
lambda = lambda0:0.1:lambda1;
k0=1./lambda;
k = k0(1):(k0(end)-k0(1))/(size(wabs,2)-1):k0(end);
wabs = interp1(k0, wabs, k);
fftnum = 4096;
res = 1/(k(1)-k(end))/2000;
res = res*size(wabs,2)/fftnum;
dr = 40;
depth = ((0:dr+2) -dr/2-1)* res;
i = FWHM / (lambda1-lambda0) / 2.355;
i = 0.5/i;
lsp = gausswin(size(wabs,2), i)';
h_sp = figure();
figure(h_sp);subplot(2,1,1);plot(lambda, lsp/max(lsp),'k', 'LineWidth', 4);
%psf of light source
lpsf = fftshift(abs(fft(lsp,fftnum)));
[v, j]=max(lpsf);
j1 = interp1(lpsf(j-3:j), j-3:j, v/2);j2 = interp1(lpsf(j:j+3), j:j+3, v/2);
fwhm1 = (j2-j1) * res;
lpsf = 20*log10(lpsf(fftnum/2-dr/2:fftnum/2+dr/2+2)/max(lpsf));    
% figure(h_sp);subplot(2,1,2);plot(depth, lpsf,'k', 'LineWidth', 2);
%cross spectrum wo reshaping
ssp = wabs.*lsp;
figure(h_sp); subplot(2,1,1);hold on;plot(lambda,ssp/max(ssp), '-.k', 'LineWidth', 4);
csp = sqrt(ssp .* lsp); 
figure(h_sp); subplot(2,1,1);hold on;plot(lambda,csp/max(csp), '--k', 'LineWidth', 4); 
%psf of cross spectrum wo reshaping
cpsf = fftshift(abs(fft(csp,fftnum)));
[v, j]=max(cpsf);
j1 = interp1(cpsf(j-3:j), j-3:j, v/2);j2 = interp1(cpsf(j:j+3), j:j+3, v/2);
fwhm2 = (j2-j1) * res;
cpsf = 20*log10(cpsf(fftnum/2-dr/2:fftnum/2+dr/2+2)/max(cpsf));    
figure(h_sp); subplot(2,1,2);hold on;plot(depth, cpsf, '--k', 'LineWidth', 2);    
%reshaped cross spectrum
csp = gausswin(size(wabs, 2), 2.2)';
figure(h_sp); subplot(2,1,1);hold on;plot(lambda,csp/max(csp),'.k', 'LineWidth', 4);
%psf of cross spectrum with reshaping
cpsf = fftshift(abs(fft(csp,fftnum)));
[v, j]=max(cpsf);
j1 = interp1(cpsf(j-3:j), j-3:j, v/2);j2 = interp1(cpsf(j:j+3), j:j+3, v/2);
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
figure(h_sp); subplot(2,1,1);hold on;plot(lambda,rsp/max(rsp),'k', 'LineWidth', 4);
%set figure titles:
subplot(2,1,1); title('spectral shape');
h_x= xlabel('wavelength(nm)');h_y=ylabel('power(a.u.)');
set(h_x, 'FontSize', 15);set(h_y, 'FontSize', 15);
% h_lg= legend('source spectrum', ...
% 'sample spectrum', 'unreshaped spectrum', ...
% 'reshaped spectrum', ...
% 'reference spectrum');
set(h_lg, 'FontSize', 15);
subplot(2,1,2); title('psf');
h_x = xlabel('depth(um)');h_y =ylabel('intensity');
set(h_x, 'FontSize', 15);set(h_y, 'FontSize', 15);
h_lg = legend(['unreshaped psf, FWHM:', num2str(fwhm2), 'um'], ...
['reshaped psf, FWHM:', num2str(fwhm3), 'um']);
set(h_lg, 'FontSize', 15);