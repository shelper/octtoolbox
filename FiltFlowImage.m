function OutputFlow=FiltFlowImage(InputFlow, FilterSize, FilterType)
%% FlowImage filtering;
FlowMask=double((InputFlow~=0));
xNoise=numel(InputFlow)/sum(FlowMask(:));
ComplexFlow=exp(1i*InputFlow);
switch FilterType
    case 'Filter: med'
        OutputFlow=medfilt2(real(ComplexFlow),FilterSize)+1i*medfilt2(imag(ComplexFlow),FilterSize);
        return; %need further consideration of phase wrapping, anyway, not frequently used.
    case 'Filter: wiener'
        [foo,Noise]=wiener2(ComplexFlow,FilterSize); 
        Noise=xNoise*Noise;%TBF?-wiener filter with variance optimized
        ComplexFlow=wiener2(ComplexFlow,FilterSize,Noise);
        OutputFlow=angle(ComplexFlow);   
    case 'Filter: gauss' 
        Filter=gausswin(FilterSize(1),2)*gausswin(FilterSize(2),2)';Filter=Filter/sum(Filter(:));
%         CorrGain=imfilter(FlowMask,Filter)+(1-FlowMask);
        OutputFlow=angle(imfilter(ComplexFlow,Filter));   
    case 'Filter: disk'
        Filter=fspecial('disk',max(FilterSize)/2);
        if size(FilterSize,2)==2
            Filter=imresize(Filter,[FilterSize(1),FilterSize(2)]);
            Filter=Filter*1/sum(Filter(:));
        end
%         CorrGain=imfilter(FlowMask,Filter)+ (1-FlowMask);
        OutputFlow=angle(imfilter(ComplexFlow,Filter));   
    case 'Filter: mean'
        Filter=fspecial('average',FilterSize);
%         CorrGain=imfilter(FlowMask,Filter)+ (1-FlowMask);
        OutputFlow=angle(imfilter(ComplexFlow,Filter));   
    otherwise
        disp('Unknown filter type.');
end

OutputFlow=OutputFlow.*FlowMask;


% SinImage=sin(InputImage*pi);
% CosImage=cos(InputImage*pi);
% switch FilterType
%     case 'Filter: med'
%         SinImage=medfilt2(SinImage,FilterSize);
%         CosImage=medfilt2(CosImage,FilterSize);            
%     case 'Filter: wiener'
%         [foo,Noise]=wiener2(SinImage,FilterSize);   
%         Noise=xNoise*Noise;%TBF?-wiener filter with variance optimized
%         SinImage=wiener2(SinImage,FilterSize,Noise);           
%         [foo,Noise]=wiener2(CosImage,FilterSize);
%         Noise=xNoise*Noise;%TBF?-wiener filter with variance optimized
%         CosImage=wiener2(CosImage,FilterSize,Noise);         
%     case 'Filter: gauss' 
%         Filter=gausswin(FilterSize(1),2)*gausswin(FilterSize(2),2)';
%         Filter=Filter/sum(Filter(:));
% %         Filter=fspecial('gaussian',FilterSize,mean(FilterSize)/3);   
%         SinImage=imfilter(SinImage,Filter);
%         CosImage=imfilter(CosImage,Filter);
%     case 'Filter: disk'
%         Filter=fspecial('disk',max(FilterSize)/2);
%         if size(FilterSize,2)==2
%             Filter=imresize(Filter,[FilterSize(1),FilterSize(2)]);
%             Filter=Filter*1/sum(Filter(:));
%         end
%         SinImage=imfilter(SinImage,Filter);
%         CosImage=imfilter(CosImage,Filter);
%     case 'Filter: mean'
%         Filter=fspecial('average',FilterSize);
%         SinImage=imfilter(SinImage,Filter);
%         CosImage=imfilter(CosImage,Filter);
%     otherwise
%         disp('Unknown filter type.');
% end
% OutputImage=angle(CosImage+1i*SinImage)/pi;    
