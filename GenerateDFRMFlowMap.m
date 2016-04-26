function AngioFlow=GenerateDFRMAgo(CaliedData, ZeroPadNum, PhaseShift)
%generate angiographic flow imaging using DPRM

    ShiftedDataPos=RampPhase(CaliedData,PhaseShift);  
    ShiftedDataNeg=RampPhase(CaliedData,-PhaseShift); 
    PosHilbertImage = abs(ifft(hilbert(ShiftedDataPos')',ZeroPadNum));
    NegHilbertImage = abs(ifft(hilbert(ShiftedDataNeg')',ZeroPadNum));

    %% use weightedratio map and intensity map
%     PosHilbertImageRatio = PosHilbertImage(2:end/2,:)./PosHilbertImage(end:-1:end/2+2,:);
%     NegHilbertImageRatio = NegHilbertImage(end:-1:end/2+2,:)./NegHilbertImage(2:end/2,:);              
%     DFRMFlowMap=PosHilbertImageRatio+NegHilbertImageRatio; 
    AngioFlow=PosHilbertImage(2:end/2,:)+NegHilbertImage(end:-1:end/2+2,:); 
    
%     DFRMFlowMap=FiltImage(DFRMFlowMap,FilterSize,FilterType);

              