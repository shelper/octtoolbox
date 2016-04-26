close all;
arr_size  = 100;
phase = randn(arr_size)*2*pi;
phase = imfilter(phase, fspecial('gaussian', 6,2) );
beam_profile = fspecial('gaussian', 100, 25);
wf = exp(1i*phase).* beam_profile;
% figure;imshow(angle(wf), []);colormap jet;
% focus = fftshift(abs(fft2(wf)));
% figure;imshow(fftshift(focus),[0, 2000]);colormap hot;

%% error reduction iteration
phase_mod = zeros(arr_size);
for m = 1:arr_size
    for n = 1:arr_size
        phase_mod(m, n) = (m+n) * pi/2;
    end
end
wf_mod = wf .* exp(1i*(phase_mod));
focus_i = abs(fft2(wf_mod));
% figure;imshow(focus_i,[]);
% foo_phase = angle(fft2(wf_mod));
focus_p = zeros(arr_size);
% focus_p(10:90, 10:90) =foo_phase(10:90, 10:90);
for n = 1:100
    wf_recon = ifft2(focus_i .* exp(1i*focus_p));
    wf_recon = exp(1i*angle(wf_recon)).* beam_profile;
    focus = fft2(wf_recon);
    focus(end/2:end, end/2:end) = 0;
    focus_p = angle(focus);
    focus_i = abs(focus);
end
figure;imshow(angle(wf_recon.* exp(-1i*(phase_mod))), []);colormap jet;
figure;imshow(angle(wf_mod.* exp(-1i*(phase_mod))), []);colormap jet;

%% add phase plate
wf1 = wf;wf2 = wf;
wf1(1:50, :)= wf1(1:50, :) * exp(1i*(pi/2));
figure;imshow(angle(wf1), []);colormap jet;
wf2(1:50, :)= wf2(1:50, :) * exp(1i*(-pi/2));
figure;imshow(angle(wf2), []);colormap jet;
focus1 = fftshift(abs(fft2(wf1)));
figure;imshow(fftshift(focus1),[0, 2000]);colormap hot;
focus2 = fftshift(abs(fft2(wf2)));
figure;imshow(fftshift(focus2),[0, 2000]);colormap hot;
focus_recon = focus1 - 1i *focus2;
wf_recon = ifft2(focus_recon);
figure;imshow(angle(wf_recon), []);colormap jet;

%% add phase modulation
phase_mod1 = zeros(arr_size);
phase_mod2 = zeros(arr_size);
for m = 1:arr_size
    for n = 1:arr_size
        phase_mod1(m, n) = (m+n) * pi/2;
        phase_mod2(m, n) = ((arr_size * 2)-(m+n)) * pi/2;
    end
end
wf1 = exp(1i*(phase+phase_mod1));
figure;imshow(angle(wf1), []);colormap jet;
wf2 = exp(1i*(phase+phase_mod2));
figure;imshow(angle(wf2), []);colormap jet;
focus1 = fftshift(abs(fft2(wf1)));
figure;imshow(fftshift(focus1),[0, 2000]);colormap hot;
focus2 = fftshift(abs(fft2(wf2)));
figure;imshow(fftshift(focus2),[0, 2000]);colormap hot;




