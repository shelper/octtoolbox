CosSpectrum=repmat(cos(1:400)',[1,400]);SinSpectrum=repmat(sin(1:400)',[1,400]);
LPM=repmat([(-50:50).^2-2500,2500-(49:-1:-49).^2]*0.25*pi,[1,2]);% figure;plot(diff(LPM))
LPM=repmat(LPM,[400,1]);
SimData=CosSpectrum.*cos(LPM)-SinSpectrum.*sin(LPM);
% figure;imshow(SimData);

fftData=ifft(SimData);
% fftData(201,:)=0;
Phase = zeros(size(SimData));
Phase(:,2:end) = angle(conj(fftData(:,2:end)).*(fftData(:,1:end-1)));

% Phase(:,1)=Phase(:,2)-(Phase(:,3)-Phase(:,2));
for i=0:3
    if mod(i,2)
        Phase(:,i*100+1:i*100+100)=-Phase(:,i*100+1:i*100+100);
    end
    Phase(:,i*100+1)=Phase(:,i*100+2)-(Phase(:,i*100+3)-Phase(:,i*100+2));
    
    if i>=1
        PhaseShift=2*Phase(:,i*100)-Phase(:,i*100-1)-Phase(:,i*100+1);
        Phase(:,i*100+1:i*100+100)=repmat(PhaseShift,[1,100])+Phase(:,i*100+1:i*100+100);
    end
end

% for i=1:3
%     if mod(i,2)
%         Phase(:,i*100+1:i*100+100)=repmat(Phase(:,i*100)-Phase(:,i*100+1)+0.5*pi,[1,100])+Phase(:,i*100+1:i*100+100);
%     else
%         Phase(:,i*100+1:i*100+100)=repmat(Phase(:,i*100)+Phase(:,i*100+1)+0.5*pi,[1,100])-Phase(:,i*100+1:i*100+100);
%     end
% end

NewfftData=abs(fftData).*exp(1i*Phase);
NewData=fft(NewfftData);
% figure;imshow(abs(NewData));
fftData=fft(NewData,[],2);
fftData(:,end/2+1:end)=0;
NewData2=ifft(fftData,[],2);
figure;imshow(20*log10(abs(fft(NewData2))/100),[-50,3]);impixelinfo;

% fftData=fft(CosSpectrum,[],2);
% fftData(:,end/2+1:end)=0;
% NewData2=ifft(CosSpectrum,[],2);
% figure;imshow(20*log10(abs(fft(NewData2))/100),[-50,3]);impixelinfo;


for i=0:3
    fftData=fft(NewData(:,i*100+1:i*100+100),[],2);
    fftData(:,[1,end/2+1:end])=0;
    NewData2(:,i*100+1:i*100+100)=ifft(fftData,[],2);
end

figure;imshow(20*log10(abs(fft(NewData2))/100),[-50,3]);impixelinfo;
figure;imshow(20*log10(abs(fft(CosSpectrum))/100),[-50,3]);impixelinfo;





