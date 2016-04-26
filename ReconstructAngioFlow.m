function CaliedData=ReconstructAngioFlow(CaliedData,PrePar, AgoPar,  ImgPar,FilePath,FileName, Procedure)
global hAngioFlowImage AngioFlow FilteredAngioFlow ;

debug=0;
% use procedure number Procedure to identify which controller is triggered, thus save processing time.    
%% reconstruction of flow with BMA correction
DeConjugate=~strcmp(PrePar.DeConjugate,'None');
if Procedure<=Const('Recon')
%     DopplerPar.Method='2ndOrder';
    switch AgoPar.Method
        case 'Speckle'
            AngioFlow=GenerateSpekAgo(CaliedData,ImgPar.ZeroPadNum);
        case 'DFRM'
            AngioFlow=GenerateDFRMAgo(CaliedData, ImgPar.ZeroPadNum,0.1*pi);
        case 'LIMP'
            AngioFlow=GenerateLIMPAgo(CaliedData, ImgPar.ZeroPadNum,0.05);
        otherwise
            errordlg('Unknown reconstruction method');
    end
end

%% Filtering and resize
if Procedure<=Const('Filter')
    FilteredAngioFlow=FiltImage(AngioFlow, ImgPar.ImgFilterSize, ImgPar.ImgFilterType);    
    FilteredAngioFlow=imresize(FilteredAngioFlow,ImgPar.ImgResize.*size(FilteredAngioFlow),'box');
end

FilteredAngioFlow=20*log10(FilteredAngioFlow);
DisplayRange=GetDisplayRange(FilteredAngioFlow, [0.5 0]);

%% Display and save
if Procedure<=Const('Display')    
    if ishandle(hAngioFlowImage)
        set(hAngioFlowImage,'CData',FilteredAngioFlow); 
    else
        figure('Name','AgoFlow','NumberTitle','off');
        hAngioFlowImage=imshow(FilteredAngioFlow,'DisplayRange',DisplayRange);impixelinfo;
    end
    if AgoPar.AutoSaveFlt
        fid = fopen([FilePath,'AgoImage',FileName(end-6:end-4),AgoPar.Method,'.dat'],'w');
        fwrite(fid, FilteredAngioFlow(21:905,:), 'float');fclose(fid);
    end 
end




