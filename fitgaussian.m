%Gaussian fitting
%FittedCurve=A*exp(-(x-xmean)^2/sigma^2)+b;
%FWHM=2*sqrt(ln2)*sigma

function estimates = fitgaussian(t, data,start_point)
% Call fminsearch with a random starting point.
%start_point = rand(1, 4);

% start_point = [0,0,0,0];

estimates = fminsearch(@expfun, start_point);
figure;plot(data);
A = estimates(1);
xmean = estimates(2);
sigma=estimates(3);
b=estimates(4);
FittedCurve = A .* exp(-(t-xmean).^2/sigma^2)+b;
hold on;plot(FittedCurve,'r')

% expfun accepts curve parameters as inputs and outputs sse,
% the sum of squares error for A * exp(-lambda * t) - Data.
    function sse = expfun(params)
        A = params(1);
        xmean = params(2);
        sigma=params(3);
        b=params(4);
        
        FittedCurve = A .* exp(-(t-xmean).^2/sigma^2)+b;
        ErrorVector = FittedCurve - data;
        sse = sum(ErrorVector .^ 2);
    end
end
