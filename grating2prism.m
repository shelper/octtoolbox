%% inputs of the system parameters
N = 1024;
alpha = 60/180*pi; % apex angle of the prism
% prism_angle = 1/180*pi;
lambda_0 = 1; lambda_n = 1.1;
d = 1000/1435;
% BK7 coefficients
% nc=[1.0396121E+00 6.0006987E-03 2.3179234E-01 2.0017914E-02 1.0104694E+00 1.0356065E+02];
% F2 coefficients
nc=[1.3453336E+00 9.9774387E-03 2.0907318E-01 4.7045077E-02 9.3735716E-01 1.1188676E+02];
% Fused Silica coefficients
% nc=[6.69422575E-01 4.34583937E-01 8.71694723E-01 4.48011239E-03 1.32847049E-02 9.53414824E+01];

%% calculate K and lambda and refractive distribution
linear_k0 = 10^6*1/lambda_n;linear_kn = 10^6*1/lambda_0;
linear_k = linear_k0:(linear_kn-linear_k0)/(N-1):linear_kn;
linear_lambda = 10^6./linear_k;lambda_c = linear_lambda(N/2);
n=sqrt(1+nc(1)*linear_lambda.^2./(linear_lambda.^2-nc(2))+nc(3)*linear_lambda.^2./(linear_lambda.^2-nc(4))+nc(5)*linear_lambda.^2./(linear_lambda.^2-nc(6)));

%% 
in_angle = asin(lambda_c/2/d);
in_angle*180/pi
beta0 = asin(linear_lambda/d - sin(in_angle));

prism_angle = 9*pi/180;
linearality = zeros(size(prism_angle, 2),N-1);
for i = 1:size(prism_angle,2)
    beta1 = beta0+prism_angle(i);
    beta1(end/2)*180/pi
    beta2 = asin(sin(beta1)./n);
    beta2(end/2)*180/pi
    beta3 = pi/2 - beta2;
    beta3(end/2)*180/pi
    beta4 = pi - alpha - beta3;
    beta4(end/2)*180/pi
    beta5 = pi/2 - beta4;
    beta5(end/2)*180/pi
    beta6 = asin(sin(beta5).*n);
    beta6(end/2)*180/pi
    beta6 = beta6 - beta6(N/2);
    beta6(end/2)*180/pi
    linearality(i,:) = diff((beta6));
end
% figure;plot(beta6)
figure;plot(linearality');
figure;plot(sin(beta6))
hold on;plot(sin(fliplr(beta0-beta0(N/2))));
figure;plot(std(linearality'));
    
    




