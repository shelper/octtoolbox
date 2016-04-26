function RawData=RemoveRef(RawData, Method)
% RawData=float(RawData);
switch Method
    case 'mean'
        RawData=RawData - repmat(mean(RawData,2),1,size(RawData,2));
    case 'median'
        RawData=RawData - repmat(median(RawData,2),1,size(RawData,2));
    case 'compmedian'
        ComplexData=fft(RawData);
        ComplexRef=median(real(ComplexData),2)+1i*median(imag(ComplexData),2);
        RawData=RawData-repmat(ifft(ComplexRef),1,size(RawData,2));
%         RawData=ifft(ComplexData);
    case 'none'
        return;
end

