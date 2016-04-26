function EvaluatePhaseStability(CaliedData,Depth)
Image=fft(CaliedData);
Image=Image(Depth,:);
PhaseLine=angle(Image(2:end).*conj(Image(1:end-1)));
figure;plot(PhaseLine,'LineStyle','none','Marker','o','MarkerEdgeColor','b');

Image=fft(CaliedData);
Image=Image(Depth*2-1,:);
PhaseLine=angle(Image(2:end).*conj(Image(1:end-1)));
hold on;plot(PhaseLine,'LineStyle','none','Marker','o','MarkerEdgeColor','r');








% Depth=48;
% Data=GetFreqComp(CaliedData,Depth*2+5,Depth*2+5);
% Phase=unwrap(angle(hilbert(Data)));
% 
% % Data1=GetFreqComp(CaliedData,1184,1204);
% % Data2=GetFreqComp(CaliedData,1228,1248);
% % Phase=unwrap(angle(hilbert(Data2)))-unwrap(angle(hilbert(Data1)));
% % figure;imshow(angle(exp(1i*diff(Phase,[],2))),[]);colormap jet;
% 
% NewData=sparse(diag(barthannwin(size(Phase,1))))*cos(Phase);
% Image=ifft(NewData);
% PhaseLine = angle(Image(Depth,2:end).*conj(Image(Depth,1:end-1)));
% figure;plot(PhaseLine);
