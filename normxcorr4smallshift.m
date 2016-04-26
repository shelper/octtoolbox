function CC=normxcorr4smallshift(InputImage1,InputImage2,ShiftX, ShiftY)
% InputImage1 is the template
ACImage1=InputImage1.^2;
ACImage2=InputImage2.^2;
AC=sqrt(sum(ACImage1(:))*sum(ACImage2(:)));

CC=zeros(ShiftX+1,ShiftY+1);
for i=-ShiftX:ShiftX
    for j=-ShiftY:ShiftY
        CCImage=InputImage1.*circshift(InputImage2,[i, j]);
        CC(i+ShiftX+1,j+ShiftY+1)=sum(CCImage(:))/AC;
    end
end
        