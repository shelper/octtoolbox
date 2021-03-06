function CaliedData=RemoveBMA(CaliedData, FullRange)
FilterSize=6;   FilterType='average';   Filter=fspecial(FilterType,FilterSize);        
% for i=1:1
%preprocessing to remove BMA
    ComplexImage=ifft(CaliedData);
    if ~strcmpi(FullRange, 'none')
        ComplexImage=fftshift(ComplexImage,1);
    else
        ComplexImage=ComplexImage(50:end/2,:);
    end
    Image=abs(ComplexImage(:,1:end-1)+ComplexImage(:,2:end))/2;
    Image=imfilter(Image,Filter);
    IntensityMask=(Image > (0.1*median(max(Image))));    
    DiffMask1=imfilter(abs(diff(ComplexImage,1,2)),Filter)./Image; 
    
    FlowComplex=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));  
    FlowComplex=exp(1i*angle(FlowComplex)).*IntensityMask;  
    
    BulkyFlow=zeros(1,size(CaliedData,2));        
    BulkyFlow(2:end)=angle(mean(FlowComplex)*size(IntensityMask,1) ...
                     ./(sum(IntensityMask)+(sum(IntensityMask)==0)));
    BulkyFlow=BulkyFlow*triu(ones(size(CaliedData,2)));
%     CaliedData=CaliedData*sparse(diag(exp(-1i*BulkyFlow)));
    
    BulkyFlow=repmat(BulkyFlow, size(CaliedData,1),1);
    CaliedData=CaliedData.*cos(BulkyFlow)-imag(hilbert(CaliedData)).*sin(BulkyFlow);

            
    ComplexImage=ifft(CaliedData);
    if ~strcmpi(FullRange, 'none')
        ComplexImage=fftshift(ComplexImage,1);
    else
        ComplexImage=ComplexImage(50:end/2,:);
    end
    DiffMask2=imfilter(abs(diff(ComplexImage,1,2)),Filter)./Image;     
    FlowMask=(DiffMask1<0.6) | (DiffMask2<0.6) & IntensityMask;% figure;imshow(FlowMask);
    FlowComplex=ComplexImage(:,2:end).*conj(ComplexImage(:,1:end-1));  
    FlowComplex=exp(1i*angle(FlowComplex)).*FlowMask;  
%     FlowComplex=FlowComplex.*FlowMask;  

    BulkyFlow=zeros(1,size(CaliedData,2));        
    BulkyFlow(2:end)=angle(mean(FlowComplex)*size(FlowMask,1) ...
                     ./(sum(FlowMask)+(sum(FlowMask)==0)));
    BulkyFlow=BulkyFlow*triu(ones(size(CaliedData,2)));
%     CaliedData=CaliedData.*(exp(-1i*BulkyFlow));
    
    BulkyFlow=repmat(BulkyFlow, size(CaliedData,1),1);
    CaliedData=CaliedData.*cos(BulkyFlow)-imag(hilbert(CaliedData)).*sin(BulkyFlow);

%% Phase Variance for flow detection
%         DiffFlow=angle(ComplexFlow(1:end-1,:).*conj(ComplexFlow(2:end,:)));    
%         VarFlow = imfilter(DiffFlow.^2,Filter)- (imfilter(DiffFlow,Filter)).^2;figure;imshow(VarFlow);
%         VarFlow = imfilter(VarFlow,Filter);
%         MaskFlow=VarFlow<0.33; 
% %         MaskFlow=circshift(MaskFlow,-FilterSize);                
%         MaskFlow=imopen(MaskFlow,strel('disk',2));
%         MaskFlow = imdilate(MaskFlow,strel('line',10,0));
        

%         MaskFlow=(VarFlow>0.3);imshow(MaskFlow);
%         FilterSize=[1,10];   FilterType='average';   Filter=fspecial(FilterType,FilterSize);    
%         VarFlow= imfilter(VarFlow,Filter);%figure;imshow((VarFlow<0.8),[]);

%         BWImage=abs(ComplexImage(1:end-1,1:end-1)+(ComplexImage(2:end,2:end)));
%         BWImage=(BWImage > (0.2*median(max(BWImage))));
%         MaskFlow= MaskFlow&BWImage;%figure;imshow(MaskFlow);
%         MaskFlow= imclose(MaskFlow,strel('disk',FilterSize));

