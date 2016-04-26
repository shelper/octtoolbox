function PhaseVar=GetWeightedPhaseVar(ComplexImage)
IntensityImage=abs(ComplexImage);
% IntensityImage=IntensityImage.*(IntensityImage>1000);
Phase=abs(angle(ComplexImage));
% MeanPhase=mean(abs(angle(sum(ComplexImage))));
MeanPhase=1.5;
% MeanPhase=MeanPhase*repmat([1,-1],[size(ComplexImage,1),(size(ComplexImage,2)+1)/2]);
PhaseVar=sum((Phase-MeanPhase).^2.*IntensityImage)./(sum(IntensityImage));
