function DispCoeff = GenDispCoeff(RawSpectrum,FilterWin)
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
DispCoeff=unwrap(angle(hilbert(SpectrumWODC)));  
foo = (DispCoeff(end)-DispCoeff(1))/(size(DispCoeff,1)-1);
DispCoeff = DispCoeff - (DispCoeff(1):foo:DispCoeff(end))';


















