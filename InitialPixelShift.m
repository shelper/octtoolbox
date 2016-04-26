function CorrectedData = InitialPixelShift(CaliedData,MaxShift)

i=1;
for p = -MaxShift : MaxShift
    for q = -MaxShift : MaxShift
        for t = -MaxShift : MaxShift
            if abs(q-p)<=MaxShift && abs(t-q)<=MaxShift && abs(t-p)<=MaxShift 
                ShiftArray(i,:)=[p, q, t];
                i = i+1;
            end
        end
    end
end

Suppression= zeros(1, size(ShiftArray,1));
% Suppression2= zeros(1, size(ShiftArray,1));
Line1=CaliedData(:,1);
for i = 1:size(ShiftArray,1)
    Line2=circshift(CaliedData(:,2),ShiftArray(i,1));
    Line3=circshift(CaliedData(:,3),ShiftArray(i,2));
    Line4=circshift(CaliedData(:,4),ShiftArray(i,3));
    Line = fft(Line2).*fft(Line3).*conj(fft(Line1).*fft(Line4));
    Line = sqrt(sqrt(Line(end/2:end)));     
%     Line = FiltImage(Line,[5,1],'Filter: mean');
    Phase=abs(angle(Line));      
    LineInt=abs(Line);   
    Suppression(i) = sum(Phase.*LineInt)/sum(LineInt);
%     Suppression1(i) = abs(angle(sum(Line)));

    Line = fft(Line1).*conj(fft(Line2));
    Line = sqrt(Line(end/2:end));
    Phase=abs(angle(Line*conj(sum(Line))));      
    LineInt=abs(Line);       
    Suppression(i) = Suppression(i) + sum(Phase.*LineInt)/sum(LineInt);
    
    Line = fft(Line2).*conj(fft(Line3));
    Line = sqrt(Line(end/2:end));
    Phase=abs(angle(Line*conj(sum(Line))));      
    LineInt=abs(Line);       
    Suppression(i) = Suppression(i) + sum(Phase.*LineInt)/sum(LineInt);    
       
    Line = fft(Line3).*conj(fft(Line4));
    Line = sqrt(Line(end/2:end));
    Phase=abs(angle(Line*conj(sum(Line))));      
    LineInt=abs(Line);       
    Suppression(i) = Suppression(i) + sum(Phase.*LineInt)/sum(LineInt);   
    
%     Line1 = fft(Line1).*conj(fft(Line2)).*conj(fft(Line1+Line4));
%     Line2 = sqrt(Line(end/2:end));     
%     Line3 = FiltImage(Line,[5,1],'Filter: mean');    
%     Phase=abs(angle(Line));      
%     LineInt=abs(Line);  
%     Suppression2(i) = sum(Phase.*LineInt)/sum(LineInt);
%     Suppression2(i) = abs(angle(sum(Line)));
%     Suppression2(i) = sum(abs(Line1+Line4-(Line2+Line3)));
end
[~, i]= min(Suppression);
CorrectedData(:,1)=Line1;
CorrectedData(:,2)=circshift(CaliedData(:,2),ShiftArray(i,1));
CorrectedData(:,3)=circshift(CaliedData(:,3),ShiftArray(i,2));
CorrectedData(:,4)=circshift(CaliedData(:,4),ShiftArray(i,3));
    
 
% figure;plot(Suppression1)
% figure;plot(Suppression2)
% figure;plot(Suppression2+Suppression1)
