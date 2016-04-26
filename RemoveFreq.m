function FilteredData=RemoveFreq(Data,Filter, Side)

fftData=fft(Data);
fftData(Filter,:)=0;
if Side==2
    fftData(end+2-Filter,:)=0;
end
FilteredData=ifft(fftData);

