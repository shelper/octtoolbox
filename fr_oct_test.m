%% system parameters
pix_num = 1024;
aline_num = 1000;
sp_profile = gausswin(pix_num);
depth_profile = [zeros(1,200),ones(1,424), zeros(1,400)]';
xpos_profile = [zeros(1,200),ones(1,400), zeros(1,400)];
oversamp_rate = 5;

%% get cplx_data;
funcs = oct_sim_funcs();
cplx_data = funcs.gen_cplx_data(pix_num, aline_num, sp_profile, depth_profile, xpos_profile,oversamp_rate);
%% add flow data
flow_data = funcs.add_flow(cplx_data, 50, -pi, 200, 400);
%% add phase alternation, const mod and sinosoidal mod
mean_phase = pi/8;
phase_amp = sin(pi/2-1.2:2.4/(pix_num-1):pi/2+1.2)';
phase_amp = phase_amp / mean(phase_amp)* mean_phase;
const_mod_data = funcs.modulate_phase(flow_data, mean_phase, 0.5);
sino_mod_data = funcs.modulate_phase(flow_data, phase_amp, 0.5);

%% 2nd order FROCT const mod at pi/8
acq_data = real(const_mod_data);
fr_data = funcs.phase_alt_deconj(acq_data, mean_phase, 0.25,4);
figure('Name','2nd order FROCT const mod at pi/8');
imshow(sqrt(fftshift(abs(fft(fr_data)), 1)),[-10, 20]);
% figure;imshow(angle(fftshift(fft(fr_data), 1)),[]);colormap jet;
% h_curve=figure();hold on;plot(mean(20*log10(fftshift(abs(fft(fr_data(:, 300:400))),1)),2),'r');

%% 2nd order FROCT sino mod at pi/8, w/o mod correction
acq_data = real(sino_mod_data);
fr_data = funcs.phase_alt_deconj(acq_data, phase_amp, 0.25, 3);
figure('Name','2nd order FROCT sino mod at pi/8, w/o mod correction');
imshow(fftshift(abs(fft(fr_data)), 1),[-10, 20]);
% figure;imshow(angle(fftshift(fft(fr_data), 1)),[]);colormap jet;
figure(h_curve);hold on;plot(mean(20*log10(fftshift(abs(fft(fr_data(:, 300:400))),1)),2),'g');

%% 2nd order FROCT sino mod at pi/8, with mod correction
acq_data = real(sino_mod_data);
fr_data = funcs.phase_alt_deconj(acq_data, phase_amp, 0.25, 2);
figure('Name','2nd order FROCT sino mod at pi/8, with mod correction');
imshow(sqrt(fftshift(abs(fft(fr_data)), 1)),[-10, 20]);
figure;imshow(angle(fftshift(fft(fr_data), 1)),[]);colormap jet;
figure(h_curve);hold on;plot(mean(20*log10(fftshift(abs(fft(fr_data(:, 300:400))),1)),2),'b');

%% reconstruct image using different methods
sino_mod_data = funcs.modulate_phase(flow_data,  phase_amp, 0.5);
sino_mod_data = funcs.modulate_phase(flow_data,  -phase_amp, 0.5);
acq_data = real(sino_mod_data);
fr_data = funcs.phase_alt_deconj(acq_data, -phase_amp, 0.25, 1);
% figure('Name','1st order FROCT by BES for +/- 0.25pi alternation');
cplx_img = fftshift(fft(fr_data),1);
figure;imshow(abs(cplx_img),[]);
figure;imshow(angle(cplx_img),[]);colormap jet;
% figure(11);hold on;plot(mean(20*log10(fftshift(abs(fft(fr_data(:, 390:410))),1)),2));
% figure;imshow(diff(angle(fftshift((fft(fr_data)), 1)),1, 2),[]);colormap jet;

fr_data = funcs.phase_alt_deconj(acq_data, phase_amp, 0.25, 2);
figure(11);hold on;plot(mean(20*log10(fftshift(abs(fft(fr_data(:, 300:400))),1)),2),'r');
figure('Name','2nd order FROCT by BES for +/- 0.25pi alternation');
imshow(fftshift(abs(fft(fr_data)), 1),[]);
% figure(11);hold on;plot(mean(20*log10(fftshift(abs(fft(fr_data(:, 300:400))),1)),2),'b');

% figure(1);imshow(angle(fftshift(fft(fr_data), 1)),[]);colormap jet;
% figure(1);imshow(fftshift(phase_diff, 1),[]);colormap jet;

fr_data = funcs.linear_mod_deconj(acq_data, phase_amp, 1);
figure('Name','1st order FROCT by LIN for +/- 0.25pi alternation');
imshow(fftshift(abs(fft(fr_data)), 1),[]);

fr_data = funcs.linear_mod_deconj(acq_data, phase_amp, 2);
figure('Name','2nd order FROCT by LIN for +/- 0.25pi alternation');
imshow(fftshift(abs(fft(fr_data)), 1),[]);

data0 = real(funcs.modulate_phase(flow_data, pi/4, 0.5));
data1 = real(funcs.modulate_phase(flow_data, pi/8, 0.5));
data2 = funcs.update_alt_phase_amp(data1, phase_amp, pi/4);

foo =(abs(fft(data0')));
figure;plot(foo(2,:))
prof0=sum((abs(fft(data0')))');
% figure;plot(prof0)
prof1=sum((abs(fft(data1')))');
prof2=sum((abs(fft(data2')))');

figure;plot([prof0', prof1', prof2']);
figure;plot(fftshift(prof1./prof0))
f


x = 1:500;
y0 = prof0(251:750);
y1 = prof1(251:750);
y2 = prof2(251:750);
% figure;plot(y)

fun0 = @(paras) sum((y0-(paras(1) * exp(-(x-251).^2/(paras(2)^2)) + paras(3))).^2)
fun1 = @(paras) sum((y1-(paras(1) * exp(-(x-251).^2/(paras(2)^2)) + paras(3))).^2)
fun2 = @(paras) sum((y2-(paras(1) * exp(-(x-251).^2/(paras(2)^2)) + paras(3))).^2)
paras0 = fminsearch(fun0,[max(y), aline_num/10, min(y0)])
paras1 = fminsearch(fun1,[max(y), aline_num/10, min(y1)])
paras2 = fminsearch(fun2,[max(y), aline_num/10, min(y2)])

gh0 =  exp(-(x-251).^2/(paras0(2)^2)) ;
gh1 =  exp(-(x-251).^2/(paras1(2)^2)) ;
gh2 =  exp(-(x-251).^2/(paras2(2)^2)) ;

figure;
% hold on;plot(prof0(251:750)
hold on;plot([gh0])
% hold on;plot(prof1(251:750))
hold on;plot([gh1])
% hold on;plot(prof2(251:750))
hold on;plot([gh2])

figure;plot(fftshift(prof1./prof2))

















figure;plot(abs(fft(funcs.update_alt_phase_amp(acq_data, phase_amp, pi/4)')))
