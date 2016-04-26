function OCTData = GetOCTData(FileName, SystemID, AlineNum, PixelNum)
debug = 0;

%% general read data program
fid = fopen(FileName);
[path, ~, ~] = fileparts(FileName);

if fid ~=-1
    switch SystemID;
        case 'datfile'
            fid = fopen(FileName, 'rb');
            OCTData = fread(fid, [PixelNum AlineNum], 'single');
        case 'TiffRaw'
            OCTData=double(imread(FileName))';
            OCTData=OCTData-65536*(OCTData>32768);
        case 'Image'
            OCTData = fread(fid, [PixelNum , AlineNum], 'single');
            OCTData=imread(FileName)';
            OCTData=double(OCTData(:,1:1000));
            fclose(fid);
        case 'THQDOCT'
            PixelNum = 2048;
            OCTData = fread(fid, [PixelNum , AlineNum], 'uint16');
%             RescalingData = load('new_sdoct.cam');   
%             RescalingData = RescalingData.phase;
            fid2=fopen([path,'/Rescale.clb'],'rb');
            RescalingData = fread(fid2, [PixelNum,1], 'float');
            fclose(fid2);
            OCTData=interp1(linspace(0, 1, PixelNum),OCTData,RescalingData, 'spline');        case 'sdOCTCam'
        case 'OCT-2000' 
            PixelNum = fread(fid, 1, 'int');   
            AlineNum = fread(fid, 1, 'int');
            fread(fid, 1, 'int');
            OCTData = fread(fid, [PixelNum, AlineNum], 'uint16');   
            RescalingData = fread(fid, [PixelNum, 1], 'single');
            RescalingData = polyval(polyfit(linspace(0, 1, PixelNum)',RescalingData,3), linspace(0, 1, PixelNum)');
            Reference = fread(fid, [PixelNum, 1], 'double');
            OCTData = bsxfun(@minus, OCTData, Reference);
            OCTData=interp1( linspace(0, 1, PixelNum)', OCTData, RescalingData, 'spline');
        case {'ssOCTCam'}
            OCTData = fread(fid, [ AlineNum,PixelNum ], 'uint16');
%             OCTData=OCTData(4:end-1,:);
%             OCTData=OCTData-2048;            
        case 'ssOCTDAQ'
            fseek(fid, 0, 'eof');FileSize = ftell(fid)/2;fseek(fid, 0, 'bof'); % get file size
            AlineNum=floor(FileSize/PixelNum);
            OCTData = fread(fid, [PixelNum, AlineNum], 'double');            
    end
    fclose(fid);
else
    OCTData = 0;
end
