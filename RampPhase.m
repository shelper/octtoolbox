function RampedData = RampPhase(CaliedData,PhaseShift) 
% digital phase ramping
if size(PhaseShift,2)==1
    PhaseShift=-PhaseShift*(1:size(CaliedData,2));
else
    PhaseShift=-PhaseShift*triu(ones(size(PhaseShift,2)),1);
end
CosSignal=sparse(diag(cos(PhaseShift*pi)));
SinSignal=sparse(diag(sin(PhaseShift*pi)));
RampedData=CaliedData* CosSignal-imag(hilbert(CaliedData))*SinSignal;
