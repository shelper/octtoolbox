p1=repmat([-1, 1],[1, 500]);
p1=p1*0.25*pi;
p2=0.7*pi*(1:1000);
p0=-0.1*pi*(1:1000);


sp=exp(1i*(p0+p1));
sn=exp(-1i*(p0+p1));
ss=sp+sn;
figure;plot(abs(fft(ss)));

scos=GetFreqComp(ss',[1,250], 2)';
ss=ss.*repmat([-1, 1],[1, 500]);
ssin=GetFreqComp(ss',[1,250], 2)';
figure;plot(real(fft(scos+1i*ssin)));
hold on;plot(imag(fft(scos+1i*ssin)),'r');


% figure;plot(real(fft(ss)));hold on;plot(imag(fft(ss)),'r');

sp=exp(1i*(p0+p1+p2));
sn=exp(-1i*(p0+p1+p2));
sd=sp+sn;
figure;plot(abs(fft(sd)));

scos=GetFreqComp(ss',[1,250], 2)';
sd=sd.*repmat([-1, 1],[1, 500]);
ssin=GetFreqComp(sd',[1,250], 2)';
figure;plot(abs(fft(scos+1i*ssin)));



sp=exp(1i*(p0+p2));
sn=exp(-1i*(p0+p2));
sd=sp+sn;
figure;plot(abs(fft(sd)));

scos=GetFreqComp(ss',[1,250], 2)';
sd=sd.*repmat([-1, 1],[1, 500]);
ssin=GetFreqComp(sd',[1,250], 2)';
figure;plot(abs(fft(scos+1i*ssin)));



figure;plot(real(fft(sd)));hold on;plot(imag(fft(sd)),'r');




fspecial('gaussian',[1,8],2)


p0=filtfilt(ones(1,8)/8,1,p0);


sp=cos(rand(phi);

