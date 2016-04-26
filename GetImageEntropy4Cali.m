function Entropy = GetImageEntropy4Cali(RawData, CaliCoeff)
    CaliK=(-1:2/2047:1)';
    OrigK=polyval([CaliCoeff, 1-CaliCoeff(1),-CaliCoeff(2)],CaliK);
    CaliedData=interp1( OrigK,RawData, CaliK, 'linear','extrap');
    Image = abs(ifft(CaliedData));
%     Entropy = sum(abs(Image(:)));
    Image = Image(1:end/2,:);
    Image=Image/sum(Image(:));
    Entropy = Image.*log(Image);
    Entropy = -sum(Entropy(:))
end