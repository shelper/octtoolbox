function PhaseError=GetFSPhaseError(AvgHRatio,InstantPhase,K,FSDriverPar)
PhaseShift=FSDriverPar(1)* sin(InstantPhase+FSDriverPar(2)) .*K;

HRatio=tan(PhaseShift)-sin(InstantPhase+FSDriverPar(2)-pi/2)*FSDriverPar(3);

PhaseError=sum((AvgHRatio-HRatio).^2);

% PhaseShift=FSDriverPar(1)* sin(InstantPhase+FSDriverPar(2)) .*K+(InstantPhase+FSDriverPar(2))*FSDriverPar(3);
hold on;plot(HRatio,'r');
