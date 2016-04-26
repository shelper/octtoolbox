function FilteredData=GetFreqComp(Data,Filter, Side, Thresh)

if nargin == 2
    Side = 2;
end
    
fftDatatmp=fft(Data);
if nargin == 4
    fftDatatmp = fftDatatmp .* (abs(fftDatatmp)>Thresh);
end
if numel(Filter)==2
    fftData=zeros(size(Data));
    if Filter(1)==1
        fftData(Filter(1),:)=fftDatatmp(Filter(1),:);
        Filter(1)=2;
    %     if iStartEnd(2)<=2
    %         FilteredData=ifft(fftData);
    %         return;
    %     end
    end
    fftData(Filter(1):Filter(2),:)=fftDatatmp(Filter(1):Filter(2),:);

    if Side==2
        fftData(end+2-Filter(2):end+2-Filter(1),:)=conj(fftDatatmp(Filter(2):-1:Filter(1),:));
    elseif Side==3
        fftData(end+2-Filter(2):end+2-Filter(1),:)=fftDatatmp(end+2-Filter(2):end+2-Filter(1),:);
    end
    FilteredData=ifft(fftData);
    return;
else
    Filter=repmat(Filter,1,size(Data,2));
    FilteredData=ifft(fftDatatmp.*Filter);
end

