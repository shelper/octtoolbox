function CaliedData=DeConjugateByFlt(CaliedData, FSPhaseShift)
% global HRatioSum;
% figure;plot(imag((fft(CaliedData,[],2)))');

AlineNum=size(CaliedData,2);
% CaliedData(:,2:2:end)=CaliedData(:,2:2:end) .* exp(1i*pi);

FSPhaseShift=repmat(-FSPhaseShift-pi,[1,AlineNum/2]);
CaliedData(:,2:2:end)=CaliedData(:,2:2:end) .* exp(1i*FSPhaseShift);
BandWidth=AlineNum/8;
CaliedData= conj(GetFreqComp(CaliedData',[1,BandWidth],3)');
