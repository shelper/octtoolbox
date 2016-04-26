function CaliedData=GenerateOCTData(SpectrumWidth, AlineNum, Position, Radius, CenterFlowRate, NoiseLevel)
CaliedData=zeros(SpectrumWidth, AlineNum);
FlowSize=Radius*2+1;
FlowPhase=zeros(FlowSize);



for i=1:FlowSize
    for j=1:FlowSize
        FlowPhase(i,j)=pi-((i - Radius-1)^2 + (j-Radius-1)^2) * pi / Radius^2;
    end
end
FlowPhase=FlowPhase.*(FlowPhase>0);
FlowPhase=FlowPhase*(triu(ones(size(FlowPhase))));

% figure;imshow(FlowPhase,[]);
% construct the OCT data of the flow flow
    I_k=zeros(1000,200);
    Iaom_k=zeros(1000,200);
    for aline=1:200
        for depth=1:500
    %                 if ((depth - flowdepth)^2 + (aline-flowaline)^2) < r^2 
                I_k(:,aline)=I_k(:,aline)+(cos(2*pi*(depth/1000:depth/1000:depth)+...
                centerspeed*sum(phase(depth,1:aline))+phasenoise(depth, aline)))';
        end
    end
    
    recon_phase= angle(fft(I_k(:,2:200)).*conj(fft(I_k(:,1:199))));
    recon_phase= recon_phase + (recon_phase>pi)*(-2*pi) +(recon_phase<-pi)*(2*pi);
    flow_p_p=recon_phase(flowdepth-90:flowdepth+90,flowaline-90:flowaline+90 );
    flow_p_p=imfilter(flow_p_p,filter);
    
%     RampRange=min(centerspeed*1.5, 0.5);

    FlowSpeed1=QuantDFRM(I_k', 1, 1, 0.02*centerspeed/0.8, 2:500);
    FlowSpeed1=imfilter(FlowSpeed1,filter);


    if loop_i==8
        figure;
        subplot(2,4,1);imshow(centerspeed*phase(flowdepth-90:flowdepth+90,flowaline-90:flowaline+90) /pi);
        subplot(2,4,5);plot(centerspeed*phase(100,10:190)/pi);
        
        subplot(2,4,2);imshow(flow_p_p/pi);colormap hot;
        subplot(2,4,6);plot(flow_p_p(90,:)/pi);
        
        subplot(2,4,3);imshow(FlowSpeed1(flowdepth-90:flowdepth+90,10:190));colormap hot;
        subplot(2,4,7);plot(FlowSpeed1(90,10:190));
        FlowSpeed2=QuantDFRM(I_k', 1, 1, 0.01*centerspeed/0.8, 2:500);
        FlowSpeed2=imfilter(FlowSpeed2,filter);
        subplot(2,4,4);imshow(FlowSpeed2(flowdepth-90:flowdepth+90,10:190));colormap hot;
        subplot(2,4,8);plot(FlowSpeed2(90,10:190));       
    end
    
%     figure;
    s0=centerspeed*phase(100,10:190)/pi;
    s1=flow_p_p(90,:)/pi-s0;
%     subplot(2,1,1);plot(s1);
    s2=FlowSpeed1(90,10:190)-s0;
%     subplot(2,1,2);plot(s2);
%     s3=FlowSpeed2(90,10:190)-s0;subplot(3,1,3);plot(s3);
    meanerror(1,loop_i)=mean(abs(s1));
    meanerror(2,loop_i)=mean(abs(s2));
    NoiseStd(1,loop_i)= std(s1);
    NoiseStd(2,loop_i)= std(s2);
    SNR(1,loop_i)=centerspeed/NoiseStd(1,loop_i);
    SNR(2,loop_i)=centerspeed/NoiseStd(2,loop_i);

%     SNR(3,loop_i)=centerspeed/NoiseStd(3,loop_i);

end

