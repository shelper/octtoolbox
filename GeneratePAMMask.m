function PAMMask=GeneratePAMMask(CaliedData,FilterSize, FilterType,ZeroPadNum)
    if IntPar.IntThresh
        CroppedIntensityImage=CroppedIntensityImage/mean(max(CroppedIntensityImage));
        IntensityMask=im2bw(CroppedIntensityImage,IntPar.IntThresh);
    else
        IntensityMask=1;
    end
end
end