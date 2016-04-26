function Var=GetWeightedVar(Phase, Intensity, Threshold)
if nargin==2
    MeanPhase=0;
else
    Mask=(Intensity>Threshold);
    MeanPhase=sum(Phase.*Mask)./(sum(Mask)+sum(Mask==0));
    MeanPhase= repmat(MeanPhase,size(Phase,1),1);
end
Var=sum((Phase-MeanPhase).^2.*Intensity)./(sum(Intensity));

