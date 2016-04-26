%% inputs of the system parameters
N = 1024;
alpha = 60/180*pi; % apex angle of the prism
% prism_angle = 1/180*pi;
lambda_0 = 1; lambda_n = 1.1; lambda_c = 1.05;
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
grating_angle = asin(lambda_c/2/d);

prism_angle = (60:70)*pi/180;
linearality = zeros(size(prism_angle, 2),N-1);
for i = 1:size(prism_angle,2)
    beta1 = asin(sin(prism_angle(i))./n);
    beta2 = pi/2 - beta1;
    beta3 = pi - alpha - beta2;
    beta4 = pi/2 - beta3;
    beta5 = asin(sin(beta4).*n);
    delta_angle = grating_angle - beta5(N/2);
    beta6 = asin(linear_lambda/d - sin(beta5+delta_angle));
    beta6 = beta6 - beta6(N/2);
    linearality(i,:) = diff(sin(beta6));
end
figure;plot(linearality');
figure;plot(std(linearality'));
    
    




