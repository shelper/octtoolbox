pixelnum=500;
alinenum=600;
FreqComp=0:0.01:1;
p=1:pixelnum;
a=1/2;
b=1.5;
phi0=0.5*pi;

RawData=zeros(pixelnum,alinenum);
RawData2=zeros(pixelnum,alinenum);
for n=1:alinenum
    phi=phi0+p/pixelnum*pi/2;
%     phi=1.6345
%     RawData(:,n)=cos(p*2*pi/10-b*sin(phi));
    RawData(:,n)=exp(1i*(p*2*pi/20+b*sin(phi0)))+exp(-1i*(p*2*pi/20+b*sin(phi0)))+exp(1i*(-p*2*pi/50+b*sin(phi0)))+exp(-1i*(-p*2*pi/50+b*sin(phi0)));
%     RawData2(:,n)=imag(hilbert(RawData(:,n)));
    phi0=phi0+a*2*pi;
end
figure;plot(abs((fft(RawData,[],2)))');
figure;plot(abs(fft(RawData)));


NewData=RawData.*repmat(sin(2*pi*a*(0:(alinenum-1))+0.5*pi),[size(RawData,1),1]);%figure;plot(abs(fft(NewData,[],2))')
H1=GetFreqComp(NewData',[1,30])';
H2=GetFreqComp(RawData',[1,30])';
% HRatio=mean(abs(H1),2)./mean(abs(H2),2);
% HRatio=sqrt(sum(H1.^2,2)./sum(H2.^2,2));
% figure;plot(HRatio);

NewData=sqrt(sum((H1(:)).^2)/sum((H2(:)).^2))*H2-1i*H1;
figure;plot((abs(fft(NewData))));%
figure;imshow(abs(fft(NewData)))

% figure;imshow(RawData,[]);

% for n=1:pixelnum 
%     RawData2(n,:)=imag(hilbert(RawData(n,:)));
% end
% RealPart=fft(RawData,[],2);
% ImagPart=fft(RawData2,[],2);
% 
% figure;plot(RealPart);
% figure;imshow(RawData)
% figure;plot(RawData(103,:));

% NewData=real((fft(RawData,[],2)))+1i*imag((fft(RawData2,[],2)));
% NewData=ifft(NewData,[],2);
% NewData=GetFreqComp(NewData',[298,306])';
% figure;imshow(abs(fft(NewData)),[]);
% figure;plot(abs(fft(NewData)));


NewData=real((fft(RawData,[],2)))+1i*real((fft(RawData2,[],2)));
NewData=ifft(NewData,[],2);
NewData=GetFreqComp(NewData',[298,306])';
figure;imshow(abs(fft(NewData)),[]);
figure;plot(abs(fft(NewData)));


RawData=zeros(pixelnum,alinenum);
for n=1:alinenum
    phi=phi0+p/pixelnum*pi/2;
%     phi=1.6345
%     RawData(:,n)=cos(p*2*pi/10+b*sin(phi));
    RawData(:,n)=(exp(1i*(p*2*pi/100+b*sin(phi0))));
    % RawData(:,n)=cos(p*2*pi/10+b*sin(phi0)+0*n*2*pi);
    phi0=phi0+a*2*pi;
end
% figure;imshow(RawData)
% figure;plot(RawData(103,:));
figure;plot(imag((fft(RawData(1:3,:),[],2)))');
hold on;plot(real((fft(RawData(1:3,:),[],2)))');





NewData=RawData.*repmat(sin(2*pi*a*(0:(alinenum-1))+0.5*pi),[size(RawData,1),1]);%figure;plot(abs(fft(NewData,[],2))')
H1=GetFreqComp(NewData',[1,50])';
NewData=RawData;
H2=GetFreqComp(NewData',[1,50])';
% NewData=sqrt(sum(abs(H1(:)))/sum(abs(H2(:))))*H2-1i*H1;
NewData=1.55*H2-1i*H1;
figure;plot((abs(fft(NewData))/1000));%imshow(abs(fft(NewData))/1000)
% figure;plot((abs(fft(NewData))/1000));%imshow(abs(fft(NewData))/1000)
figure;imshow(abs(fft(NewData))/1000)

NewData=RawData.*repmat(sin(2*pi*a*(0:(alinenum-1)))+sin(6*pi*a*(0:(alinenum-1))),[size(RawData,1),1]);%figure;plot(abs(fft(NewData,[],2))')
H1=GetFreqComp(NewData',1,30)';
NewData=RawData.*repmat(1+sin(4*pi*a*(0:(alinenum-1))+0.5*pi)+sin(8*pi*a*(0:(alinenum-1))+0.5*pi),[size(RawData,1),1]);
H2=GetFreqComp(NewData',1,30)';
NewData=sqrt(sum((H1(:)).^2)/sum((H2(:)).^2))*H2-1i*H1;
figure;plot((abs(fft(NewData))/1000));%imshow(abs(fft(NewData))/1000)
NewData=sqrt(sum(abs(H1(:)))/sum(abs(H2(:))))*H2-1i*H1;
figure;plot((abs(fft(NewData))/1000));%imshow(abs(fft(NewData))/1000)
figure;imshow(abs(fft(NewData))/1000)