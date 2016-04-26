 %from oct 2000
color = ['r','g','b','y','k'];
max_psf1 = 0;max_psf2 = 0;
N = 5;

for dn = 1:5
    depth = dn *100;    
    k0 = (0:1/1023:1);
    ks = k0.^2 * 0.156 + k0*0.844 ;
%     disp = 2*pi* 500*k0 * 0.156 - 2*pi* 500*k0.^2*0.156;
    disp = 0;
    sig = cos(2*pi* depth *ks+disp);
    
    foo_cs = rand(size(sig))>0.5;
    sum(foo_cs)
    ind1 = 1:1024; ind2 = 1:1024;
    ind1 = ind1(foo_cs==0);ind2 = ind2(foo_cs==1);
%     sig(ind1) = interp1(ind2, sig(ind2), ind1);
    sig = sig(ind2);
    ks = ks(ind2);

%     foo_cs = [zeros(1, 500), ones(1, length(sig)-500)];
%     sig = sig .* foo_cs ;
    sig_calied = spline(ks,sig, k0);
%     rsn = 2048;
%     k0 = (0:1/(rsn-1):1);
% %     hold on;plot(sig,'r');
%     [win, ind, psf_corr] = CalibrateK(ks, k0, N, 'besselwin');
%     wn = size(win, 1);
%     sig_calied = zeros(1, length(k0));
%     for i = 1 : length(k0)
%         sig_calied(i) = sum(win(:, i)'.* sig(ind(i):ind(i)+wn-1));
%     end
%     sig_calied = sig_calied .* gausswin(rsn)';
%     sig_linear = cos(2*pi* depth *k0) .* gausswin(rsn)';
%     
% %     disp = spline(ks, disp, k0);
%     sig_calied = sig_calied .*exp(-1i * disp);
%     
    psf = 20*log10(abs(fft(sig_calied)));
    psf = psf(1:512);
%     psf = psf + 20*log10(psf_corr);
    max_psf1 = max(max_psf1,max(psf));
%     psf = angle(fft(sig_calied));
%     psf = psf +2*pi * (psf<0);
    figure(1);hold on;plot(psf-max_psf1,['-.',color(mod(dn,5)+1)],'linewidth',3);
    
    psf = 20*log10(abs(fft(sig_linear)));
    psf = psf(1:512);
    max_psf2 = max(max_psf2,max(psf));
%     psf = angle(fft(sig_linear));
%     psf = psf +2*pi * (psf<0);
    figure(1);hold on;plot(psf-max_psf2,['',color(mod(dn,5)+1)],'linewidth',3)  
end

ylim([-100, 0]);
xlim([0, 512]);

% r = length(k0)/length(ks);
% t = pi/(r*(r-0.5)) * N/2;
% g = exp(-(-5:5).^2 * r^2/(2*t));
% figure;plot(g)
% gain = abs(fft(g, rsn));



