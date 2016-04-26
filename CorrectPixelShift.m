function CaliedData=CorrectPixelShift(CaliedData,MaxShift)
% debug = 1;

[PixelNum,AlineNum]=size(CaliedData);
%% decide the sign of the begin phase method 2
% Line1=CaliedData(:,end);
% Line2=CaliedData(:,end-1);
% % ShiftSteps=-MaxShift:MaxShift;
% SuppressionPos=zeros(MaxShift*2+1,1);
% SuppressionNeg=zeros(MaxShift*2+1,1);
% for m=-MaxShift:MaxShift
%     LinePos=Line2+1i*(circshift(Line1,m) ./ sin(PhaseShift) - Line2./tan(PhaseShift));
%     LineNeg=Line2-1i*(circshift(Line1,m) ./ sin(PhaseShift) - Line2./tan(PhaseShift));
%     LinePos=log(abs(fft(LinePos)));
%     LineNeg=log(abs(fft(LineNeg)));
%     SuppressionPos(m+MaxShift+1)=sum(abs(LinePos));
%     SuppressionNeg(m+MaxShift+1)=sum(abs(LineNeg));
% end
% [foo, Index1]=min(SuppressionPos(1:end));  
% [foo, Index2]=min(SuppressionNeg(1:end)); 
% 
% if SuppressionPos(Index1) > SuppressionNeg(Index2) 
%     PhaseShift=-PhaseShift;
%     TotalShift= Index2 - MaxShift - 1;
% else
%     TotalShift= Index1 - MaxShift - 1;
% end
%% decide the sign of the begin phase method 1
% % if PhaseShift
% foo=fft(CaliedData(:,1 : 100)); %get the complex image by fourier transform of RawData, select 100 Alines in the middle to minimize the computation 
% foo = foo(2:end/2,:); %select one of the conjugate half 
% PhaseDiff=angle(foo (:,2:2:end).*conj(foo (:,1:2:end-1))); %calculate the phase difference between Odd Alines and Even Alines (Even Phase – Odd Phase)
% BinaryImage= abs(foo(:,2:2:end)) +abs(foo(:,1:2:end-1));
% BinaryImage=(BinaryImage>500);
% % get BinaryImage as a mask for PhaseDiff based on the Intensity:
% PhaseDiff = PhaseDiff .*BinaryImage;% screen the PhaseDiff ro remove low intensity phase noise
% PhaseDiff =sum(PhaseDiff (:)); %calculate the averaged phase difference of above
% if PhaseDiff > 0
%    CaliedData = circshift(CaliedData,[0 -1]);
% %     AvgPhase = - AvgPhase; %change the sign of HRatio if the phase difference is negative            
% end
% PhaseShift= - atan(FRHRatio)*2;
%% phase jittering correction by glass cover
% ACSignalDepth=568;
% ACSignal=fft(CaliedData);
% ACSignal=ACSignal(ACSignalDepth-5: ACSignalDepth+5,:);
% PixelShift=zeros(AlineNum,1);
% PhaseShiftPerPixel=ACSignalDepth*2*pi/PixelNum;
% RefSignal=conj(ACSignal(:,1));
% for n=1:AlineNum  
%     PhaseShift=angle(sum(ACSignal(:,n).*RefSignal));
% 	PixleShiftArray=[PhaseShift-2*pi, PhaseShift, PhaseShift+2*pi]/PhaseShiftPerPixel;
%     JitError=abs(PixleShiftArray-round(PixleShiftArray));
%     [foo,i]=min(JitError);
%     PixelShift(n)=round(PixleShiftArray(i));   
%     CaliedData(:,n)=circshift(CaliedData(:,n),PixelShift(n));
% end
% % define the start phase signal sign
% foo = fft(CaliedData(:,AlineNum/2-49 : AlineNum/2+50)); %get the complex image by fourier transform of RawData, select 100 Alines in the middle to minimize the computation 
% foo = foo(2:end/2,:); %select one of the conjugate half 
% PhaseDiff=foo (:,2:2:end).*conj(foo (:,1:2:end-1));%calculate the phase difference between Odd Alines and Even Alines (Even Phase – Odd Phase)
% PhaseDiff=angle(sum(PhaseDiff(:)));
% if PhaseDiff < 0
%    CaliedData = circshift(CaliedData,[0 1]);
% end
%% correct by minimizing diff between Alines
% ACSignal=GetFreqComp(CaliedData,ACSignalDepth,2);
% AvgACSignal=ACSignal(IndexRange,1);
% 
% for n=2:AlineNum
%     AlineSpectrum=ACSignal(:,n);    
%     JitError=1000000000;    
%     for m = -MaxShift:MaxShift
%         if JitError > sum(abs(AvgACSignal-AlineSpectrum(IndexRange+m)))
%             JitError = sum(abs(AvgACSignal-AlineSpectrum(IndexRange+m)));
%             PixelShift = m;
%         end
%     end    
%     AvgACSignal = AvgACSignal*(n-1)/n + AlineSpectrum(IndexRange+PixelShift)/n; 
%     CaliedData(IndexRange,n)=CaliedData(IndexRange+ PixelShift,n);
% end

% if Step == 2
%     JitterNoise=1000000000000;
%     Line1=sum(CaliedData(:,1:Step:end),2);
%     Line2=sum(CaliedData(:,2:Step:end),2);
%     for m=-MaxShift:MaxShift  
%         Error=sum(abs(Line1-circshift(Line2,m)));
%         if JitterNoise > Error
%             JitterNoise = Error;
%             Shift=m;
%         end
%     end    
%     CaliedData(:,2:2:end)=circshift(CaliedData(:,2:2:end),[Shift,0]);
% end

%     
% Step=2;
% for n=AlineNum:-Step:AlineNum-500
%     Line1=ACSignal(:,n);
%     Line2=ACSignal(:,n-Step);
%     JitError=1000000000;    
%     for m=-MaxShift:MaxShift  
%         if JitError > sum(abs(Line1-circshift(Line2,m)))
%             JitError = sum(abs(Line1-circshift(Line2,m)));
%             Shift=m;
%         end
%     end    
%     ACSignal(:,n-Step)=circshift(ACSignal(:,n-Step),Shift);
% %     CaliedData(:,n-Step)=circshift(CaliedData(:,n-Step),Shift);
% end
% ACSignal=mean(ACSignal(:,AlineNum:-2:AlineNum-500),2);
% ACSignal=repmat(ACSignal,[1,AlineNum]);
% 
% JitError=zeros(MaxShift*2+1,AlineNum);
% for m=-MaxShift:MaxShift  
%     TempData=circshift(CaliedData, [m,0]); 
%     TempData=GetfreqComp(TempData,MirrorLocation,2);
% %     TempData=abs(fft(TempData))-abs(fft(ACSignal));
%     JitError(m+MaxShift+1,:)=sum(abs(TempData-ACSignal));
% %     JitError(m+MaxShift+1,:)=sum(TempData(MirrorLocation(1):MirrorLocation(2),:));    
% end
% [foo,PixShift]=min(JitError);
% PixShift=PixShift-MaxShift-1;
% PixShift=min(0,PixShift);
% 
% for i=1:AlineNum
%     CaliedData(:,i)=circshift(CaliedData(:,i),PixShift(i));
% end

% CaliedData=CaliedData-ACSignal;
%% phase jittering correction based on phase
% ShiftSteps=-MaxShift:MaxShift;
% Suppression=zeros(MaxShift*2+1,1);
% PhaseShift=-atan(FRHRatio)*2;
% % ComplexImage = fft(CaliedData);
% % ComplexImage = ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));
% % PhaseImage = angle(ComplexImage);
% % MeanPhase=0;
% for n=AlineNum:-1:2  
%     Line1=CaliedData(:,n);
%     Line2=CaliedData(:,n-1).*exp(1i*PhaseShift);
%     for m=ShiftSteps
%         Line = fft(circshift(Line2,m)).*conj(fft(Line1));
%         Line = Line(end/2:end);
%         Phase=angle(Line);
% %         Int=abs(Line(2:end/2));
% %         Suppression(m-ShiftSteps(1)+1) = sum(sum(Phase(2:end/2).*abs(Line(2:end/2)));
% %         Suppression(m-ShiftSteps(1)+1) = sum(abs(Phase).*abs(Line));
% %         Suppression(m-ShiftSteps(1)+1) = sum(abs(Phase));
% %         MeanPhase =  sum(sum(Phase.*Int)/sum(Int));
%         Suppression(m-ShiftSteps(1)+1) = sum(abs(Phase).*abs(Line));
%     end
% %     n
%     [foo, Index]=min(Suppression(1:size(ShiftSteps,2)));  
%     CaliedData(:,n-1)=circshift(Line2,ShiftSteps(Index));
%     PhaseShift=-PhaseShift;
%     n
% end
%% phase jittering correction based on intensity by Even Alines shift
% Suppression= zeros(MaxShift*2+1, AlineNum);
% ComplexData= zeros(PixelNum,AlineNum);
% CaliedData(:,1:2:end-1) = circshift(CaliedData(:,1:2:end-1),-MaxShift);
% for m= 1: MaxShift*2+1
%     if m == 1
%         CaliedData(:,1:2:end-1) = circshift(CaliedData(:,1:2:end-1),-MaxShift);
%     else
%         CaliedData(:,1:2:end-1) = circshift(CaliedData(:,1:2:end-1),1);
%     end
%     for n = AlineNum-1:-2:2
%         ComplexData(:,n+1) = CaliedData(:,n)+1i*(CaliedData(:,n+1) ./ sin(PhaseShift) - CaliedData(:,n)./tan(PhaseShift));        
%         ComplexData(:,n) = CaliedData(:,n)+1i*(CaliedData(:,n-1) ./ sin(PhaseShift) - CaliedData(:,n)./tan(PhaseShift));        
%     end
%     Suppression(m,:) = sum(log(abs(fft(ComplexData))));
% end
% CaliedData(:,1:2:end-1) = circshift(CaliedData(:,1:2:end-1),-MaxShift);
% Suppression(:,1:2:end-1)=flipud(Suppression(:,1:2:end-1));
% 
% TotalShift = 0;
% for n = AlineNum-1:-1:1
%     StartShift=max(-MaxShift, -MaxShift+TotalShift)+MaxShift+1;
%     EndShift=min(MaxShift, MaxShift+TotalShift)+MaxShift+1;    
%     [foo, Index]=min(Suppression(StartShift:EndShift, n+1));
%     TotalShift=TotalShift - MaxShift + (Index - StartShift);
%     if abs(TotalShift) <= MaxShift
%         CaliedData(:,n)= circshift(CaliedData(:,n),TotalShift);
%     else
%         TotalShift = 0; 
%     end
% end
%% phase jittering correction by complex phase    
% Line0=0;
% StartShift=-MaxShift;
% EndShift=MaxShift;
% PrevShift=0;
% StartPixel=1;
% EndPixel=1376;
% AvgPhaseShift=mean(PhaseShift(StartPixel:EndPixel));
% for n=AlineNum:-1:2
%     Line1=CaliedData(StartPixel:EndPixel,n);
%     Line2=CaliedData(StartPixel:EndPixel,n-1);
%     MinConjNoise= 100000000;
%     for m = StartShift : EndShift
%         ComplexLine=Line1+1i*(circshift(Line2,m)./sin(PhaseShift)-Line1./tan(PhaseShift) );
% %         CaliedData(:,n)+1i*(CaliedData(:,n+1) ./ sin(PhaseShift) - CaliedData(:,n)./tan(PhaseShift));        
%         Line =  (abs(fft(ComplexLine)));
%         ConjNoise=sum(Line)+sum(abs(Line-Line0));            
% %         Phase=angle(Line);
%         if MinConjNoise > ConjNoise
%            MinConjNoise = ConjNoise;
%            Shift=m;
%         end
%     end
%     CaliedData(:,n-1)=circshift(CaliedData(:,n-1),Shift);   
%     
% %     if n<AlineNum-1
% %         ComplexLine=CaliedData(StartPixel:EndPixel,n+2) ...
% %             +1i*(CaliedData(StartPixel:EndPixel,n+1)./sin(PhaseShift) ... 
% %             -CaliedData(StartPixel:EndPixel,n+2)./tan(PhaseShift) );
% %         Line0 = log (abs(fft(ComplexLine))) ;
% %     end
% %     
% %     if abs(Shift)>abs(PrevShift);
% %         StartShift=max(-MaxShift, -MaxShift+Shift);
% %         EndShift=min(MaxShift, MaxShift+Shift);
% %         PrevShift=MaxShift;
% %     end
% %     n
%     PhaseShift=-PhaseShift;    
% 
% end
%% phase jittering by log abs line
% StartShift=-MaxShift;
% EndShift=MaxShift;
% Filter = ones(5,1)/5;
% % PrevShift=0;
% % Shift = 0;
% for n=AlineNum:-1:2
%     Line1=CaliedData(:,n);
%     Line2=CaliedData(:,n-1);
%     ConjNoise=0;
%     for m=StartShift:EndShift
%         Line=Line1+1i*(circshift(Line2,m) ./ sin(PhaseShift) - Line1./tan(PhaseShift));
%         Line=(abs(fft(Line)));
%         Line = (imfilter(Line.^2, Filter)-(imfilter(Line,Filter)).^2);
%         if ConjNoise < sum(Line)
% %            StdImg=(imfilter(Line.^2, Filter)-(imfilter(Line,Filter)).^2)./imfilter(AbsImg, Filter);
%             ConjNoise = sum(Line);
%             Shift=m;
%         end
%     end
% %     TotalShift=TotalShift+Shift;
% % 
% %     if (Shift~=0 && n==AlineNum) 
% %         StartShift=max(-MaxShift, -MaxShift+Shift);
% %         EndShift=min(MaxShift, MaxShift+Shift);
% %     else
% %         StartShift=-MaxShift;
% %         EndShift=MaxShift;
% %     end  
%     
%     CaliedData(:,n-1)=circshift(Line2,Shift);
%     PhaseShift=-PhaseShift;
% %     PrevShift=Shift;
%     n
% 
% end
%% Phase Jittering correction 
% ConjNoise=zeros(1, MaxShift*2+1);   
% ShiftSteps=-MaxShift:MaxShift;
% Mask=1;
% for n=1:AlineNum-1  
%     Line1=CaliedData(:,n);
%     Line2=CaliedData(:,n+1);
%     for m=ShiftSteps
%         Line=Line1+1i*(circshift(Line2,m) ./ sin(PhaseShift) - Line1./tan(PhaseShift));
%         
%         Line=abs(fft(Line));
% %         Mask=(Line(2:end/2)>300) & (Line(end:-1:end/2+2)>300);
%         ConjNoise(m-ShiftSteps(1)+1)=sum(Line(150:end-148));
%     end
% %     [~, Index]=min(ConjNoise(1:size(ShiftSteps,0)));  
% %     CaliedData(:,n+1)=circshift(CaliedData(:,n+1),ShiftSteps(Index));
% %     Line=abs(fft(Line1+1i*(CaliedData(:,n+1) ./ sin(PhaseShift) - Line1./tan(PhaseShift))));
% %     Line=imfilter(Line, ones(10,1));
% %     Mask=(Line(2:end/2)>Line(end:-1:end/2+2))-(Line(2:end/2)<Line(end:-1:end/2+2));
% %     Mask=Mask.*((Line(2:end/2)>3000)|Line(end:-1:end/2+2)>3000);   
%     PhaseShift=-PhaseShift;
% end  
% correct phase shift by glass plate
% if ByGlassPlate
%     for n=1:AlineNum-1  
%         [~,Index]= max(xcorr(CaliedData(1233:1239,n), CaliedData(1233:1239,n+1),MaxShift, 'coeff'));
%         PixelShift= Index-(MaxShift+1);
%         CaliedData(:,n+1)= circshift(CaliedData(:,n+1),PixelShift);
%     end
% else
%     CaliedData(:,2:2:end) = CaliedData(:,2:2:end) .* repmat(exp(1i*PhaseShift),[1,size(CaliedData,2)/2]); 
% pixel shift by Cost function
% phase jittering correction by phase minimization for odd & even  
% Phase0=0;
% ComplexLine0=1;
% for n=1:AlineNum-2
%     Line1=CaliedData(:,n);
%     Line2=CaliedData(:,n+2);  
%     ComplexLine=abs(fft(Line1))+abs(fft(Line2));
%     Mask=(abs(ComplexLine(2:end/2))>(1000));
%     for m=ShiftSteps
%         ComplexLine=conj(fft(Line1)).*fft(circshift(Line2,m));
%         Phase=angle(ComplexLine(2:end/2).*ComplexLine0);
%         Phase=abs(Phase(Mask));
%         ConjNoise(m-ShiftSteps(1)+1)=sum(Phase);
%     end
%     [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));  
%     CaliedData(:,n+2)=circshift(CaliedData(:,n+2),ShiftSteps(Index));        
%     ComplexLine0=fft(Line1).*conj(fft(CaliedData(:,n+2)));
%     ComplexLine0=ComplexLine0(2:end/2);
% %     Phase0=angle(ComplexLine(2:end/2));
% end
% Line1=CaliedData(:,end/2);
% Line2=CaliedData(:,end/2+1);  
% ComplexLine=abs(fft(Line1))+abs(fft(Line2));
% Mask=(abs(ComplexLine(2:end/2))>(1000));
% for m=ShiftSteps
%     ComplexLine=conj(fft(Line1)).*fft(circshift(Line2,m));
%     Phase=angle(ComplexLine(2:end/2));
%     Phase=abs(Phase(Mask)-median(PhaseShift));
% %     figure;plot(Phase);
%     ConjNoise(m-ShiftSteps(1)+1)=sum(Phase);
% end
% [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));
% CaliedData(:,1:2:end)=circshift(CaliedData(:,1:2:end),ShiftSteps(Index));
% phase jittering correction by phase minimization
% process the first 2 Alines
% Line1=conj(fft(CaliedData(:,1)));
% Line1=Line1(2:end/2);
% Mask=(abs(Line1)>300);
% for m=ShiftSteps
%     Line2=fft(circshift(CaliedData(:,2),m));
%     Line2=Line2(2:end/2);
%     Phase=angle(Line2.*Line1);
%     Phase=Phase(Mask)-AvgPhase;
%     Phase=abs(Phase);    
%     ConjNoise(m-ShiftSteps(1)+1)=sum(Phase);
% end
% [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));  
% CaliedData(:,2)=circshift(CaliedData(:,2),ShiftSteps(Index));
% process the rest Alines with bulky motion correction
% AvgPhase=-AvgPhase;
% foo=1;
% for n=2:AlineNum-1  
%     Line1=conj(fft(CaliedData(:,n)));
%     Line1=Line1(2:end/2);
% %     Line2=fft(CaliedData(:,n+1));      
% %     Line2=Line2(2:end/2);
% %     Mask=(abs(Line1)+abs(Line2)>1000);
%     Mask=(abs(Line1)>300);    
%     for m=ShiftSteps
%         Line2=fft(circshift(CaliedData(:,n+1),m));
%         Line2=Line2(2:end/2);
%         Phase=angle(Line2.*Line1);
%         Phase=Phase(Mask);
% %         Phase=Phase(1:end-10);
%         Phase=Phase-AvgPhase;
% %         figure;plot(Phase);
%         Phase=abs(Phase); 
%         ConjNoise(m-ShiftSteps(1)+1)=sum(Phase(1:end-10));
%     end
%     [~, Index]=max(ConjNoise(1:size(ShiftSteps,2)));  
%     CaliedData(:,n+1)=circshift(CaliedData(:,n+1),ShiftSteps(Index));
%     foo=fft(CaliedData(:,n-1)).*conj(fft(CaliedData(:,n+1)));
%     foo=foo(2:end/2);
%     foo=sqrt(foo);    
% %     foo=fft(CaliedData(:,n)).*conj(fft(CaliedData(:,n+1)));
% %     foo=foo(2:end/2)*exp(1i*AvgPhase);
% %     foo=conj(foo);
%     AvgPhase=-AvgPhase;
% end
    %%test method2
    % for n=1:2
    %     Line1=conj(fft(CaliedData(:,n)));
    %     Line2=fft(CaliedData(:,n+1));      
    %     Mask=(abs(Line1(2:end/2))+abs(Line2(2:end/2))>1000);
    %     for m=ShiftSteps
    %         Line2=fft(circshift(CaliedData(:,n+1),m));
    %         ComplexLine=Line2.*Line1;
    %         Phase=angle(ComplexLine(2:end/2));
    %         Phase=Phase(Mask)-AvgPhase;
    % %         figure;plot(Phase);
    %         Phase=abs(Phase);    
    %         ConjNoise(m-ShiftSteps(1)+1)=sum(Phase);
    %     end
    %     [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));  
    %     CaliedData(:,n+1)=circshift(CaliedData(:,n+1),ShiftSteps(Index));
    %     AvgPhase=-AvgPhase;
    % end
    % for n=3:AlineNum-1  
    %     Line0=fft(CaliedData(:,n-2)) .* conj(fft(CaliedData(:,n-1)));
    %     Line0=Line0(2:end/2);       
    %     Line1=conj(fft(CaliedData(:,n)));
    %     Line1=Line1(2:end/2);
    %     Mask=(abs(Line1)>1000);    
    %     for m=ShiftSteps
    %         Line2=fft(circshift(CaliedData(:,n+1),m));
    %         Line2=Line2(2:end/2);
    %         Phase=angle(Line2.*Line1.*Line0);
    %         Phase=Phase(Mask);
    %         Phase=Phase-mean(Phase);
    %         Phase=abs(Phase); 
    %         ConjNoise(m-ShiftSteps(1)+1)=sum(Phase);
    %     end
    %     [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));  
    %     CaliedData(:,n+1)=circshift(CaliedData(:,n+1),ShiftSteps(Index));
    % end
    
    %%test method1
    % for n=1:AlineNum-1  
    %     Line1=conj(fft(CaliedData(:,n)));
    %     Line2=fft(CaliedData(:,n+1));      
    %     Mask=(abs(Line1(2:end/2))+abs(Line2(2:end/2))>1000);
    %     for m=ShiftSteps
    %         Line2=fft(circshift(CaliedData(:,n+1),m));
    % %         ComplexLine=Line2.*Line1.*Line1.*Line0;
    % %         ConjNoise(m-ShiftSteps(1)+1)=sum(abs(imag(ComplexLine)));       
    %         ComplexLine=Line2.*Line1;
    %         ComplexLine0=Line2.*Line0;
    %         Phase=angle(ComplexLine(2:end/2))-angle(ComplexLine0(2:end/2))/2;
    %         Phase=Phase(Mask)-AvgPhase;
    % %         figure;plot(Phase);
    %         Phase=abs(Phase);    
    %         ConjNoise(m-ShiftSteps(1)+1)=sum(Phase);
    %     end
    %     [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));  
    %     CaliedData(:,n+1)=circshift(CaliedData(:,n+1),ShiftSteps(Index));
    %     Line0=CaliedData(:,n);
    %     AvgPhase=-AvgPhase;
    % end
% Phase jittering correction by minimizing BG noise.
% Line0=fft(CaliedData(:,1));
% Line0=Line0(567:572);
% for n=1:AlineNum-1  
%     for m=ShiftSteps
%         Line=fft(circshift(CaliedData(:,n+1),m));   
%         Line=Line(567:572);
%         %% suppresion ratio indicator
%         ConjNoise(m-ShiftSteps(1)+1)=sum(abs(real(Line0/n-Line))+abs(imag(Line0/n-Line)));
%     end
%     [~, Index]=min(ConjNoise(1:size(ShiftSteps,2)));  
%     CaliedData(:,n+1)=circshift(CaliedData(:,n+1),ShiftSteps(Index));
%     Line0=Line0+Line;
% %     ComplexLine=conj(fft(CaliedData(:,n))).*fft(CaliedData(:,n+1));
% %     BMAPhase=angle(ComplexLine(2:end/2));
% %     Mask=(abs(ComplexLine(2:end/2))>(1000000));
% %     BMAPhase=mean(BMAPhase(Mask))-median(PhaseShift);
% %     CorrPhaseShift=BMAPhase-PhaseShift;
% %     PhaseShift=-PhaseShift;
% %         ShiftSteps=max(-MaxShift, -MaxShift + ShiftSteps(Index)): min(MaxShift, MaxShift + ShiftSteps(Index));
% end
% phase jittering correction using the conjugate suppresion ratio         
% PixelShift=zeros(1,AlineNum);
% failed to correct the phase jittering by applying conj phase shift for + and - images
%             CaliedDataN = CaliedData;
%             CaliedDataP = CaliedData;
% 
%             CaliedDataN(:,2:2:end) = CaliedDataN(:,2:2:end) .* repmat(exp(1i*PhaseShift),[1,size(CaliedData,2)/2]);
%             CaliedDataP(:,2:2:end) = CaliedDataP(:,2:2:end) .* repmat(exp(-1i*PhaseShift),[1,size(CaliedData,2)/2]);
% 
%             PixelShiftN=zeros(1,size(CaliedDataN,2)-2);
%             PixelShiftP=zeros(1,size(CaliedDataN,2)-2);
%             for i=1:size(CaliedDataN,2)-1
%                 LineN1=fft(CaliedDataN(:,i));
%                 LineN1(end/2+1)=abs(LineN1(end/2+1));      
%                 LineN1(1)=0;      
%                 LineN1(end/2:-1:2)= conj(LineN1(end/2+2:end));
% 
%                 LineN2=fft(CaliedDataN(:,i+1));
%                 LineN2(end/2+1)=abs(LineN2(end/2+1));
%                 LineN2(1)=0;      
%                 LineN2(end/2+2:end)=conj(LineN2(end/2:-1:2));
% 
%                 AvgIntensity=abs(LineN1)+abs(LineN2);
%                 Threshold=max(AvgIntensity)/50; 
%                 AvgIntensity= AvgIntensity.*(abs(LineN1)>Threshold).*(abs(LineN2)>Threshold);
% 
%                 RawLineN1= ifft(AvgIntensity.*exp(1j*angle(LineN1)));    
%                 RawLineN2= ifft(AvgIntensity.*exp(1j*angle(LineN2)));
% 
%                 [~,Index]= max(xcorr(RawLineN1, RawLineN2, MaxShift, 'coeff'));
%                 PixelShiftN(i)= Index-(MaxShift+1);        
%                 CaliedDataN(:,i+1)= circshift(CaliedDataN(:,i+1),PixelShiftN(i));
% 
% 
%                 LineP1=fft(CaliedDataP(:,i));
%                 LineP1(end/2+1)=abs(LineP1(end/2+1));      
%                 LineP1(1)=0;      
%                 LineP1(end/2+2:end)=conj(LineP1(end/2:-1:2));
% 
%                 LineP2=fft(CaliedDataP(:,i+1));
%                 LineP2(end/2+1)=abs(LineP2(end/2+1));
%                 LineP2(1)=0;      
%                 LineP2(end/2+2:end)=conj(LineP2(end/2:-1:2));
% 
%                 AvgIntensity=abs(LineP1)+abs(LineP2);
%                 Threshold=max(AvgIntensity)/50; 
%                 AvgIntensity= AvgIntensity.*(abs(LineP1)>Threshold).*(abs(LineP2)>Threshold);
%                 RawLineP1= ifft(AvgIntensity.*exp(1j*angle(LineP1)));    
%                 RawLineP2= ifft(AvgIntensity.*exp(1j*angle(LineP2)));
% 
%                 [~,Index]= max(xcorr(RawLineP1, RawLineP2, MaxShift, 'coeff'));
%                 PixelShiftP(i)= Index-(MaxShift+1);        
%                 CaliedDataP(:,i+1)= circshift(CaliedDataP(:,i+1),PixelShiftP(i));
%             end
% correct phase shift by correlation       
%     for n=1:AlineNum-1  
%         Line1=fft(CaliedData(:,n));
%         Line1(end/2+1)=abs(Line1(end/2+1));      
%         Line1(1)=0;      
%         Line1(end/2+2:end)=conj(Line1(end/2:-1:2));
%         Line2=fft(CaliedData(:,n+1));
%         Line2(end/2+1)=abs(Line2(end/2+1));
%         Line2(1)=0;      
%         Line2(end/2+2:end)=conj(Line2(end/2:-1:2));
% 
%         AvgIntensity=abs(Line1)+abs(Line2);
%         Threshold=max(AvgIntensity)/100; 
%         AvgIntensity= AvgIntensity.*(abs(Line1)>Threshold).*(abs(Line2)>Threshold);
%         RawLine1= ifft(AvgIntensity.*exp(1j*angle(Line1)));    
%         RawLine2= ifft(AvgIntensity.*exp(1j*angle(Line2)));
% 
%         [~,Index]= max(xcorr(RawLine1, RawLine2,MaxShift, 'coeff'));
%     % 	[~,Index]= max(xcorr(CaliedData(:,i), CaliedData(:,i+1),MaxShift, 'coeff'));
%         PixelShift= Index-(MaxShift+1);
%         CaliedData(:,n+1)= circshift(CaliedData(:,n+1),PixelShift);
%         TempData(:,n+1)= circshift(TempData(:,n+1),PixelShift);
%     end   
% correct phase shift by correlation, aligning A-lines by odd and even A-lines respectively
%     if ~isempty(PhaseShift)
%         foo=fft(CaliedData);
%         foo=foo(end/20:end/2-end/20,:);
%         foo=angle(foo(:,2:2:end).*conj(foo(:,1:2:end-1)));
%         foo=sum(foo(:));
%         if foo>0
%             PhaseShift=-PhaseShift;
%         end
%     end  
%     PreviousSum=0;
%     for i=1:size(CaliedData,2)-2
%         Line1=fft(CaliedData(:,i));
%         Line2=fft(CaliedData(:,i+2));
%         AvgIntensity=abs(Line1)+abs(Line2);
%         AvgIntensity= AvgIntensity.*(AvgIntensity> (max(AvgIntensity)/2));
%         RawLine1= ifft(AvgIntensity.*exp(1j*angle(Line1)));    
%         RawLine2= ifft(AvgIntensity.*exp(1j*angle(Line2)));
% 
%         [~,Index]= max(xcorr(RawLine1, RawLine2,MaxShift, 'coeff'));
%         PixelShift= Index-(MaxShift+1);
%         CaliedData(:,i+2)= circshift(CaliedData(:,i+2),PixelShift);
%         
%         %find lines where Intensity is majorly focused at high freq.
%         AvgIntensity([1:round(end/5), end+2-round(end/5):end])=0;
%         IntensitySum=sum(AvgIntensity);
%         if IntensitySum > PreviousSum
%             TempIndex=i;
%             PreviousSum=IntensitySum;
%         end     
%     end
%     if mod(TempIndex,2)
%         TempIndex=TempIndex+1;
%     end
%     CaliedData(:,TempIndex)=CaliedData(:,TempIndex) .*cos(PhaseShift)-imag(hilbert(CaliedData(:,TempIndex))).*sin(PhaseShift);
%     [~,Index]= max(xcorr(CaliedData(:,TempIndex-1), CaliedData(:,TempIndex),MaxShift, 'coeff'));
%     PixelShift= Index-(MaxShift+1);
%     CaliedData(:,TempIndex)=CaliedData(:,TempIndex) .*cos(PhaseShift)+imag(hilbert(CaliedData(:,TempIndex))).*sin(PhaseShift);    
%     CaliedData(:,2:2:end)=circshift(CaliedData(:,2:2:end),PixelShift);      

% end
% retrieve the CaliedData by the saved TempData
% CaliedData=TempData;


% for i=2:Step:size(CaliedData,2)-Step
%     Line1=fft(CaliedData(:,i));
%     Line2=fft(CaliedData(:,i+Step));
%     Phase=angle(Line2.*conj(Line1));
%     Line1=abs(Line1);
%     Line2=abs(Line2);
%     AvgIntensity=Line1+Line2;
%     Threshold=max(AvgIntensity)/50; 
%     
%     AvgIntensity= AvgIntensity.*(Line1>Threshold).*(Line2>Threshold);
%     RawLine1=ifft(AvgIntensity);    
%     RawLine2= ifft(AvgIntensity.*exp(1j*Phase));
%     
% 	[~,Index]= max(xcorr(RawLine1, RawLine2,MaxShift, 'coeff'));
%     PixelShift= Index-(MaxShift+1);
%     TempData(:,i+Step)= circshift(TempData(:,i+Step),PixelShift);
% end

% figure;plot(PixelShift);

% function CaliedData=CorrectPixelShift2(CaliedData,ByDiff) 
% %% self pixelshift correction  
% % ComplexImage=ifft(CaliedData);
% % ComplexImage=ComplexImage(2:end/2,:);
% % ComplexImage=conj(ComplexImage(:,2:end)).*ComplexImage(:,1:end-1);
% % DopplerFlow = angle(ComplexImage);
% % DopplerIntensityImage = sqrt(abs(ComplexImage));
% % 
% % DopplerFlow=DopplerFlow.*(DopplerIntensityImage>0.5)/pi;
% % DopplerFlow(end/2:end,:)=0;DopplerFlow(1:100,:)=0;
% % DopplerFlow=DopplerFlow-median(mean(DopplerFlow)*687./sum(DopplerFlow~=0))*(DopplerFlow~=0);
% % 
% % PixelShift=DopplerFlow./repmat((1:687)'/687,[1,999]);
% % % figure;plot(DopplerFlow(:,200));
% % PixelShift=round(sum(PixelShift)./sum(PixelShift~=0));
% % % figure;plot(PixelShift)
% % for i=1:size(CaliedData,2)-1
% %      CaliedData(:,i+1)= circshift(CaliedData(:,i+1),sum(PixelShift(1:i)));
% % end          
% 
% %%pixelshift correction by AC signal
% ComplexImage=fft(CaliedData);
% [~, PeakPosition]=max(sum(abs(ComplexImage(250:300,:)),2));
% PeakPosition=PeakPosition+250-1;
% % PeakPosition=210;
% % Image=abs(fft(CaliedData));  figure;imshow(Image/10000);
% 
% if ByDiff
% %         tic;
% %     ACData=zeros(size(CaliedData));
% % %         ACData(PeakPosition+1,:)=ComplexImage(PeakPosition+1,:);
% % %         ACData(size(ACData,1)-PeakPosition+1,:)=ComplexImage(size(ACData,1)-PeakPosition+1,:);
% %     ACData(PeakPosition-1:PeakPosition+1,:)=ComplexImage(PeakPosition-1:PeakPosition+1,:);
% %     ACData(size(ACData,1)-PeakPosition+1:size(ACData,1)-PeakPosition+3,:)=ComplexImage(size(ACData,1)-PeakPosition+1:size(ACData,1)-PeakPosition+3,:);
% %     ACData=ifft(ACData);
%     PeakPosition=[PeakPosition-1 PeakPosition+1];
%     ACData=GetFreqComp(CaliedData,PeakPosition);
%     DiffN3= sum(abs(ACData(3:end-4,1:end-1)-ACData(6:end-1,2:end)));
%     DiffN2= sum(abs(ACData(2:end-3,1:end-1)-ACData(4:end-1,2:end)));
%     DiffN1= sum(abs(ACData(2:end-2,1:end-1)-ACData(3:end-1,2:end)));
%     Diff0=  sum(abs(ACData(:,1:end-1)-ACData(:,2:end)));    
%     DiffP1= sum(abs(ACData(3:end-1,1:end-1)-ACData(2:end-2,2:end)));
%     DiffP2= sum(abs(ACData(4:end-1,1:end-1)-ACData(2:end-3,2:end)));
%     DiffP3= sum(abs(ACData(6:end-1,1:end-1)-ACData(3:end-4,2:end)));
%     [~,PixelShift]=min([DiffN3;DiffN2;DiffN1;Diff0;DiffP1;DiffP2;DiffP3]);
%     PixelShift=PixelShift-4;
% %         toc;
% else
% %         tic;
%     ComplexImage=fft(CaliedData);  
%     Phase=angle(ComplexImage(PeakPosition,2:end).*conj(ComplexImage(PeakPosition,1:end-1)))/pi;
%     PixelShift=round(Phase*size(CaliedData,1)/PeakPosition/2);
% %         toc;
% end
%     figure;plot(PixelShift);  
% 
% for i=1:size(CaliedData,2)-1
%      CaliedData(:,i+1)= circshift(CaliedData(:,i+1),sum(PixelShift(1:i)));
% end          
% CaliedData([1:5,end-4:end],:)=0.01;
% pixel correction by maximizing the correlation
% function CaliedData=CorrectPixelShift3(CaliedData)
% MaxShift=5;
% % Filter=zeros(size(CaliedData));
% % Range1=284;
% % Filter([Range1,end+2-Range1],:)=1;
% % fftData=fft(CaliedData);
% % FixNoiseData=ifft(fftData.*Filter);
% % FixNoiseData=sparse(diag(barthannwin(size(FixNoiseData,1))))*FixNoiseData;%Spectrum Reshaping
% 
% PixelShift=zeros(1,size(CaliedData,2));
% for i=1:size(CaliedData,2)-1
% 	[~,Index]= max(xcorr(CaliedData(:,i), CaliedData(:,i+1),MaxShift, 'coeff'));
% function CaliedData=CorrectPixelShift(CaliedData)
% 
%% pixel correlation by phase
ComplexImage=ifft(CaliedData);  
[~, PeakPosition]=max(sum(abs(ComplexImage(50:end-50,:)),2));
PeakPosition = PeakPosition+49;
ComplexImage=ComplexImage(PeakPosition,:);
Phase=angle(conj(ComplexImage(2:end)).*ComplexImage(1:end-1));
Phase = Phase * (PixelNum /2 + 1)  / (pi *(PeakPosition-1)); 
PixelShift=round(Phase);
PixelShift(abs(PixelShift) == 2) = -PixelShift(abs(PixelShift) == 2);

for i=1:size(CaliedData,2)-1
     CaliedData(:,i+1)= circshift(CaliedData(:,i+1),sum(PixelShift(1:i)));
end          

%     Phase=angle(conj(ComplexImage(2:end)).*ComplexImage(1:end-1))/pi;
%     PhaseP=Phase+2;
%     PhaseN=Phase-2;
%     PixelShift=[Phase*size(CaliedData,1)/(PeakPosition-1)/2;
%                 PhaseP*size(CaliedData,1)/(PeakPosition-1)/2;
%                 PhaseN*size(CaliedData,1)/(PeakPosition-1)/2];
%     Error=abs((PixelShift-round(PixelShift))./round(PixelShift));
%     [~,MinErrorIndex]=min(Error);
%     PixelShift=round(PixelShift);
%     PixelShift=PixelShift(MinErrorIndex+(0:998)*3);
%     figure;plot(PixelShift);
%  
%% phase jittering by odd and even line
% StartShift=-MaxShift;
% EndShift=MaxShift;
% 
% for n=AlineNum:-1:AlineNum-1
%     Line1=CaliedData(:,n);
%     Line2=CaliedData(:,n-1);
%     ConjNoise=1000000000;    
%     for m=StartShift:EndShift
%         Line=Line1+1i*circshift(Line2,m);
%         Line=(abs(fft(Line)));
%         if ConjNoise > sum(Line)
%             ConjNoise = sum(Line);
%             Shift=m;
%         end
%     end    
%     CaliedData(:,n-1)=circshift(Line2,Shift);
% end
% 
% for n=AlineNum-2:-1:2
%     Line1=CaliedData(:,n-1);
%     Line2=conj(fft(CaliedData(:,n)));
%     Line3=conj(fft(CaliedData(:,n+1)));
%     Line4=(fft(CaliedData(:,n+2)));
%     
%     ConjNoise=10000000;
%     for m=StartShift:EndShift
%         Line=fft(circshift(Line1,m)).*Line3;
%         if ConjNoise > sum(abs(angle(Line))); 
%             ConjNoise = sum(abs(angle(Line)));
%             Shift=m;
%         end
%     end
%     CaliedData(:,n-1)=circshift(Line1,Shift);
%     n
% end
%% Phase Jittering correction by phase ramp fitting
% ComplexImage = fft(CaliedData);
% ComplexImage = ComplexImage(1:PixelNum/2+1,2:end).*conj(ComplexImage(1:PixelNum/2+1,1:end-1));
% % ComplexImage = FiltImage(ComplexImage,[3, 1], 'Filter: mean');
% PhaseJittering=(0:PixelNum/2)'*[-MaxShift:MaxShift]/PixelNum*2*pi;  
% IntensityImage = sqrt(abs(ComplexImage));
% Var=zeros(MaxShift*2+1,AlineNum-1);
% for i=1:MaxShift*2+1     
%     NewComplexImage = ComplexImage .* exp(-1i* repmat(PhaseJittering(:,i), 1,AlineNum-1));
%     PhaseVar(i,:) = GetWeightedPhaseVar(NewComplexImage);
% end
% [~, iMinVar]=min(PhaseVar);
% iMinVar=iMinVar-MaxShift-1;
% % figure;plot(iVar);
% for i = 1 : AlineNum-1
%     CaliedData(:,i+1)=circshift(CaliedData(:,i+1),sum(iMinVar(1:i)));
% end

%% Phase jittering correction with doppler correction
% PhaseMod=0;
% Suppression=zeros(2*MaxShift+1,1);
% if PhaseMod
%     CaliedData(:,1:4) = InitialPixelShift(CaliedData(:,1:4),MaxShift);
%     StepNum=2*MaxShift+1;
%     Suppression=zeros(StepNum,1);
%     for n=5:AlineNum
%         Line1=CaliedData(:,n-3);
%         Line2=CaliedData(:,n-2);
%         Line3=CaliedData(:,n-1);
%         Line4=CaliedData(:,n);
%         for m=-MaxShift:MaxShift
%             Line = fft(Line2).*fft(Line3).*conj(fft(Line1).*fft(circshift(Line4,m)));
%             Line = sqrt(sqrt(Line(end/2+200:end-50)));     
%             Line = FiltImage(Line,[5,1],'Filter: mean');             
%             Phase=abs(angle(Line));      
%             LineInt=abs(Line);  
%             Suppression(m+MaxShift+1) = sum(Phase.*LineInt)/sum(LineInt);
%         end
%         [foo, Index]=min(Suppression);  
%         CaliedData(:,n)=circshift(Line4,Index-MaxShift-1);
%     end
% else
%     CaliedData(:,1:3) = NonModInitialPixelShift(CaliedData(:,1:3),MaxShift);
%     for n=4:AlineNum
%         Line1=CaliedData(:,n-2);
%         Line2=CaliedData(:,n-1);
%         Line3=CaliedData(:,n);        
%         for m=-MaxShift:MaxShift
%             Line = fft(Line2).^2.*conj(fft(Line1).*fft(circshift(Line3,m)));
%             Line = sqrt(sqrt(Line(end/2:end)));
%             %Line = FiltImage(Line,[5,1],'Filter: mean');   
%             Phase=abs(angle(Line));      
%             LineInt=abs(Line);   
%             Suppression(m+MaxShift+1) = sum(Phase.*LineInt)/sum(LineInt);
%         end
%         [~, Index]=min(Suppression);  
%         CaliedData(:,n)=circshift(Line3,Index-MaxShift-1);
%     end
% end

