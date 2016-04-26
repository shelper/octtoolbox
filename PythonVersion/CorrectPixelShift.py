function [RawData,PhaseShift]=CorrectPixelShift(RawData,MaxShift,PhaseShift,ByGlassPlate)
debug = 1;

[PixelNum,AlineNum]=size(RawData);
    
%% correct phase shift
if ByGlassPlate
    RawData=GetFreqComp(RawData,[286,289],2);
    for i=1:AlineNum-1
    	[~,Index]= max(xcorr(RawData(:,i), RawData(:,i+1),MaxShift, 'coeff'));
        PixelShift= Index-(MaxShift+1);
        RawData(:,i+1)= circshift(RawData(:,i+1),PixelShift);
    end
else
    Suppression=zeros(1, MaxShift*2+1);   
    ShiftSteps=-MaxShift:MaxShift;
    Line0=0;
    for n=1:AlineNum-1  
        Line1=RawData(:,n);
        Line2=RawData(:,n+1);      
        for m=ShiftSteps
            Line=abs(fft(Line1+1i*(circshift(Line2,m) ./ sin(PhaseShift) - Line1./tan(PhaseShift))));
            Suppression(m-ShiftSteps(1)+1)=sum(abs(Line-Line0)/2)-sum(abs(Line(2:end/2)-Line(end:-1:end/2+2)));
        end

        [~, Index]=min(Suppression(1:size(ShiftSteps,2)));        
        RawData(:,n+1)=circshift(RawData(:,n+1),ShiftSteps(Index));        
        Line0=abs(fft(Line1-1i*(circshift(Line2,ShiftSteps(Index)) ./ sin(PhaseShift) - Line1./tan(PhaseShift))));
        ShiftSteps=max(-MaxShift, -MaxShift + ShiftSteps(Index)): min(MaxShift, MaxShift + ShiftSteps(Index));
    end
end
