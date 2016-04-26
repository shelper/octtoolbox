function CombineBackForwardFrames(SystemID)
[FileNameArray,FilePath]=uigetfile({'*.*';'*.raw';'*.bin';'*.mat'}, 'Pick raw data files','MultiSelect', 'on');
if FilePath
    FileNameArray=char(FileNameArray);
    FileNum=size(FileNameArray,1);
    ProgressBar = waitbar(0,'Combine BackForward Frames...');

    for FileIndex=1:2:FileNum
        FileName=deblank(FileNameArray(FileIndex,:));
        ForwardFrame=uint16(GetOCTData([FilePath,FileName],SystemID));
        FileName=deblank(FileNameArray(FileIndex+1,:));
        BackFrame=fliplr(uint16(GetOCTData([FilePath,FileName],SystemID)));
        imwrite([ForwardFrame,BackFrame],[FilePath,FileName(1:end-4),'BF.tif']);  
        waitbar(FileIndex/FileNum,ProgressBar);    
    end
    close(ProgressBar);
end