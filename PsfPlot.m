function psf = PsfPlot(FilePath, SystemID, DispCoeff, AlineNum, PixelNum)
RawData = GetOCTData(FilePath, SystemID, AlineNum, PixelNum);

pixel_norm = (linspace(0, 1, PixelNum).^2 - linspace(0, 1, PixelNum)) * 4;
disp = exp(1i * pixel_norm * DispCoeff);
OCTData = bsxfun(@times, RawData, disp);

OCTData = bsxfun(@times, OCTData, gausswin(PixelNum));

psf = mean(abs(fft(OCTData, 2048)),2);
psf = 20*log10(psf(1:1024));
psf = psf - max(psf);

figure();plot(psf);



