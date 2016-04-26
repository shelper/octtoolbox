function OutputImage=FiltImage(InputImage, FilterSize, FilterType)
%% FlowImage filtering;    
ImageMask=double((InputImage~=0));
xNoise=numel(InputImage)/numel(InputImage(InputImage~=0));
% if bFlowImage
%     InputImage=cos(InputImage*pi)+i*sin(InputImage*pi);
% end
switch FilterType
    case 'Filter: med'
        OutputImage=medfilt2(real(InputImage),FilterSize)+1i*medfilt2(imag(InputImage),FilterSize);
    case 'Filter: wiener'
        [foo,Noise]=wiener2(InputImage,FilterSize); 
        Noise=xNoise*Noise;%TBF?-wiener filter with variance optimized
        OutputImage=wiener2(InputImage,FilterSize,Noise);
    case 'Filter: gauss' 
        Filter=gausswin(FilterSize(1),2)*gausswin(FilterSize(2),2)';
        Filter=Filter/sum(Filter(:));
%         Filter=fspecial('gaussian',FilterSize,mean(FilterSize)/3);   
%         CorrGain=imfilter(ImageMask,Filter);
%         CorrGain=CorrGain+(~ImageMask);
        OutputImage=imfilter(InputImage,Filter);%./CorrGain;
    case 'Filter: disk'
        Filter=fspecial('disk',max(FilterSize)/2);
        if size(FilterSize,2)==2
            Filter=imresize(Filter,[FilterSize(1),FilterSize(2)]);
            Filter=Filter/sum(Filter(:));
        end
% %         CorrGain=imfilter(ImageMask,Filter);
% %         CorrGain=CorrGain+(~ImageMask);
        OutputImage=imfilter(InputImage,Filter);%./CorrGain;
    case 'Filter: mean'
        Filter=fspecial('average',FilterSize);
%         CorrGain=imfilter(ImageMask,Filter);
%         CorrGain=CorrGain+(~ImageMask);
        OutputImage=imfilter(InputImage,Filter);%./CorrGain;
    otherwise
        disp('Unknown filter type.');
end
% OutputImage=OutputImage.*ImageMask;

% if bFlowImage
%     OutputImage=angle(OutputImage)/pi;
% end

