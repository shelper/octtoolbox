function [win, ind, corr] = CalibrateK(ks, k0, N, method)
    % ks -- the sampled k
    % data -- the interferogram
    % k0 -- the linear k
    % N -- the window size  
%     ks = (ks -ks(1))/(ks(end)-ks(1)); % normalize ks to [0,1]
    dk0 = 1/(length(k0)-1);
    osr = length(k0)/length(ks);
    ind = zeros(1, length(k0));
    wn = N-2;
    win = zeros(wn, length(k0));
    g = zeros(1, wn);
    for m = 1:length(k0)
        [~,fooi] = sort(abs(k0(m)-ks),'ascend');
        i = min(fooi(1:wn)); % find out the index of the first ks within window to k0(i)
        ind(m) = i;
        for n = 1:wn
            kn = (k0(m) - ks(i+n-1)) /dk0;
            g(n) = gridwin(kn, osr, N, method);
            win(n,m) = g(n);
        end
%         figure;plot(g)
    end
    corr = get_corr(osr,N,length(ks)/2, length(k0), method);
end

function g = gridwin(k,r,N, method)
    if strcmp(method,'gausswin')
        if abs(k) <= N
            t = pi/(r*(r-0.5)) * N/2;
            g = exp(-k^2*r^2/(2*t));
        else
            g = 0;
        end
    elseif strcmp(method,'besselwin')
        if abs(k) <= N/2
            beta = pi*sqrt(N^2*(r-0.5)^2/r^2 - 0.8);
            g = besseli(0,beta * sqrt(1-(2*k/N)^2))/N;
        else
            g = 0;
        end
    end
end

function psf_corr = get_corr(r,N,dmax, Ms, method)
    psf_corr = zeros(1, dmax);
    if strcmp(method,'gausswin')
        t = pi/(r*(r-0.5)) * N/2;
        for i = 1:dmax
            depth = i-1;
            psf_corr(i) = exp((depth/(2*dmax))^2*4*t)/sqrt(4*t);
        end
%         figure;plot(psf_corr)
    elseif strcmp(method,'besselwin')
        beta = pi*sqrt(N^2*(r-0.5)^2/r^2 - 0.8);
        for i = 1:dmax
            depth = i-1;
            foo = sqrt((depth*pi*N/Ms)^2 - beta^2);
            foo = foo/sin(sqrt((depth*pi*N/Ms)^2 - beta^2));
            psf_corr(i) = foo;
        end
    end
        
end


