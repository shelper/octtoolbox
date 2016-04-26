close all;
clear all;
delta_x = 0.015;%mm, x100um
aline_num = 1024;
aline_rate = 100000; %lines/s
width_psf = 0.02;%mm, x1000um
vc = delta_x*aline_rate; %mm/s
n1_level = 0.0;
n2_level = 0.0;

% am = angle(randn(1,40) + 1i * randn(1,40));
% vm = rand(1,40)*40;
am = linspace(-pi, pi, 40);
vm = ones(1, 40)*40;
for n = 1:40
    corr_n1 = randn(1, aline_num)*n1_level;
    corr_n2 = randn(1, aline_num)*n2_level;
    corr_bg =  exp(-(vc* (1/aline_rate)).^2/width_psf.^2) + corr_n1; 
    corr_bg = corr_bg .* (corr_bg<1) .* (corr_bg>0) + 0.01*((corr_bg>1) + (corr_bg<0));
    v_scan = 1i*vc*exp(-1i* linspace(0, (aline_num-1)/aline_num * 2*pi, aline_num)) + vm(n) *exp(1i*am(n));
    delta_x = abs(v_scan) * (1/aline_rate); 
    corr_x = exp(-delta_x.^2/width_psf.^2) + corr_n1 + corr_n2;
    corr_x = corr_x .* (corr_x<1) .* (corr_x>0) + 0.01*((corr_x>1) + (corr_x<0));
    diff_corr = corr_x - corr_bg;
%     diff_corr = 40*50000*((width_psf *sqrt(log(1./corr_x))).^2-(width_psf * sqrt(log(1./corr_bg))).^2) ;
%     figure;plot(diff_corr);
    A(n) = mean(diff_corr(1:aline_num/4));
    B(n) = mean(diff_corr(aline_num/4+1:aline_num/2));
    C(n) = mean(diff_corr(aline_num/2+1:aline_num*3/4));
    D(n) = mean(diff_corr(aline_num*3/4+1:aline_num));
    %% use A, B
%     am_calc(n) = atan(B(n)/A(n));
%     if A(n)<0 && B(n)<0
%         am_calc(n) = am_calc(n)-pi ;
%     elseif  A(n)<0 && B(n)>0 
%         am_calc(n) = am_calc(n)+pi ;
%     end
%     vm_calc(n) = abs(A(n)/(2*vc*sin(am_calc(n))));
%     am_calc(n) = angle(exp(1i*(am_calc(n)+pi/4)));
    %% use B, C
%     am_calc(n) = atan(B(n)/-C(n));
%     if B(n)<0 && C(n)>0
%         am_calc(n) = am_calc(n)-pi ;
%     elseif  C(n)>0 && B(n)>0 
%         am_calc(n) = am_calc(n)+pi ;
%     end
%     vm_calc(n) = abs(C(n)/(2*vc*sin(am_calc(n))));
%     am_calc(n) = angle(exp(1i*(am_calc(n)+pi/4)));
    %% use C, D
%     am_calc(n) = atan(D(n)/C(n));
%     if D(n)>0 && C(n)>0
%         am_calc(n) = am_calc(n)-pi ;
%     elseif  C(n)>0 && D(n)<0 
%         am_calc(n) = am_calc(n)+pi ;
%     end
%     vm_calc(n) = abs(C(n)/(2*vc*sin(am_calc(n))));
%     am_calc(n) = angle(exp(1i*(am_calc(n)+pi/4)));
    %% use D, A
%     am_calc(n) = atan(D(n)/-A(n));
%     if D(n)>0 && A(n)<0
%         am_calc(n) = am_calc(n)-pi ;
%     elseif  A(n)<0 && D(n)<0 
%         am_calc(n) = am_calc(n)+pi ;
%     end
%     vm_calc(n) = abs(A(n)/(2*vc*sin(am_calc(n))));
%     am_calc(n) = angle(exp(1i*(am_calc(n)+pi/4)));
    %% use A-C and B-D
    %% Get Diff Phase Variation for tracking
end
MoveVel = ((B-D) + 1i * (A-C))*1300 ;
MoveVel = MoveVel* exp(1i*3*pi/4);

figure;plot(MoveVel(1:end-3),'-+');
figure;plot(cumsum(MoveVel(1:end))/100, '-+');

% am_calc = angle(exp(1i*(am_calc+pi/4)));
% figure;plot(C,B,'-o');
% figure;plot(am_calc,'ro-');hold on;plot(am,'bo-');
% figure;plot(vm_calc,'ro-');hold on;plot(vm,'bo-');

% figure;plot(1:40, abs(vm_calc));
