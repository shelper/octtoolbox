function GenCaliCoeff(RawSpectrum,FilterWin, FileName)
RawSpectrumFFT=fft(RawSpectrum);
SpectrumWODC = zeros(size(RawSpectrumFFT));
SpectrumWODC([FilterWin(1):FilterWin(2),end+2-FilterWin(2):end+2-FilterWin(1)])= ...
            RawSpectrumFFT([FilterWin(1):FilterWin(2),end+2-FilterWin(2):end+2-FilterWin(1)]);
% RawSpectrumFFT(abs(RawSpectrumFFT)<(max(abs(RawSpectrumFFT))/100))=0;
SpectrumWODC=ifft(SpectrumWODC);
% save([FilePath,'f_',FileName], 'SpectrumWODC','-ascii')
% figure;plot(abs(RawSpectrumFFT));pause;
% FreqRange= input('pls input the frequency range: ');
% RawSpectrumFFT([1:FreqRange(1),FreqRange(2):SpectrumWidth+2-FreqRange(2),SpectrumWidth+2-FreqRange(1):end])=0;
CaliCoeff=unwrap(angle(hilbert(SpectrumWODC))); 
CaliCoeff = (CaliCoeff -CaliCoeff(1)) /(CaliCoeff(end)-CaliCoeff(1));

fid = fopen(FileName, 'wb');
fwrite(fid, CaliCoeff, 'double');
fclose(fid);

% global CaliCoeff;
% CaliCoeff=[];
% RawSpectrumFFT=fft(RawSpectrum);
% RawSpectrumFFT([1:end/10,end/2-round(end/10):end/2+2+round(end/10),end:-1:end+2-end/10])=0;
% RawSpectrumFFT(end/2-end/10:end/2+2-end/10)=0;
% % RawSpectrumFFT(abs(RawSpectrumFFT)<(max(abs(RawSpectrumFFT))/100))=0;
% SpectrumWODC=ifft(RawSpectrumFFT);
% % save([FilePath,'f_',FileName], 'SpectrumWODC','-ascii')
% % figure;plot(abs(RawSpectrumFFT));pause;
% % FreqRange= input('pls input the frequency range: ');
% % RawSpectrumFFT([1:FreqRange(1),FreqRange(2):SpectrumWidth+2-FreqRange(2),SpectrumWidth+2-FreqRange(1):end])=0;
% 
% SpectrumWidth=max(size(SpectrumWODC));
% Phase=unwrap(angle(hilbert(SpectrumWODC))); 
% Step=(Phase(SpectrumWidth)-Phase(1))/(SpectrumWidth-1);
% for i=1:SpectrumWidth
%     LinearPhase(i)=Phase(1)+Step*(i-1);
% end
% 
% LeftPointIndex(1)=1;
% LeftPointCoeff(1)=1;
% RightPointIndex(1)=2;
% RightPointCoeff(1)=0;
% LeftPointIndex(SpectrumWidth)= SpectrumWidth-1;
% LeftPointCoeff(SpectrumWidth)= 0;
% RightPointIndex(SpectrumWidth)= SpectrumWidth;
% RightPointCoeff(SpectrumWidth)= 1;
% j=2;
% for i=2:SpectrumWidth
%     while LinearPhase(i)>Phase(j)
%         j=j+1;
%     end
%     LeftPointIndex(i)=j-1;
%     RightPointIndex(i)=j;
%     LeftPointCoeff(i)=(Phase(j)-LinearPhase(i))/(Phase(j)-Phase(j-1));
%     RightPointCoeff(i)=(LinearPhase(i)-Phase(j-1))/(Phase(j)-Phase(j-1));  
% end
% 
% CaliCoeff(1,:)=LeftPointIndex;
% CaliCoeff(2,:)=LeftPointCoeff;
% CaliCoeff(3,:)=RightPointIndex;
% CaliCoeff(4,:)=RightPointCoeff;
% 
% % save([FilePath,'CaliCoeff.mat'],'CaliCoeff');
