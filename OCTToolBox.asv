function varargout = OCTToolBox(varargin)
%OCTToolBox M-file for OCTToolBox.fig
%      OCTToolBox, by itself, creates a new OCTToolBox or raises the existing
%      singleton*.
%
%      H = OCTToolBox returns the handle to a new OCTToolBox or the handle to
%      the existing singleton*.
%
%      OCTToolBox('Property','Value',...) creates a new OCTToolBox using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to OCTToolBox_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      OCTToolBox('CALLBACK') and OCTToolBox('CALLBACK',hObject,...) call the
%      local function named CALLBACK in OCTToolBox.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OCTToolBox

% Last Modified by GUIDE v2.5 29-Sep-2011 12:06:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OCTToolBox_OpeningFcn, ...
                   'gui_OutputFcn',  @OCTToolBox_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function OCTToolBox_OpeningFcn(hObject, eventdata, handles, varargin)
global CaliCoeff;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for OCTToolBox
handles.output = hObject;

% Enhance handles structure
guidata(hObject, handles);
% TBF?-load dispersion and calibration files here
CaliCoeff=[];
GetParameters(handles);
function varargout = OCTToolBox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;
%% Main Panel
function SystemID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SystemID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FileIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%% PreProcess
function ReferenceMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DispFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DispFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DeConjOrder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DeConjWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeConjWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% Intensity
function IntScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgFilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DopplerIntThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DopplerIntThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function IntDispRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ImgFilterType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgFilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ImgFilterSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgFilterSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% ImageEnhance and enhancement
function ZeroPadNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgFilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ImageFOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end%% Main Panel
function ImgResize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgResize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% Flow
function DopplerFilterSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DopplerFilterSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DopplerFilterType_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DopplerMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DopplerMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% angio
function AgoThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function AgoMethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Real Processing Functions                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Golbal functions
function SystemID_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.SystemID=get(handles.SystemID,'String');  
PrePar.SystemID=PrePar.SystemID{get(handles.SystemID,'Value')};
function LoadFiles_Callback(hObject, eventdata, handles)
global FilePath FileNameArray Ref PrePar;
[FileNameArray,FilePath]=uigetfile({'*.*';'*.raw';'*.bin';'*.mat'}, 'Pick raw data files','MultiSelect', 'on');
FileNameArray=char(FileNameArray);
FileNum=size(FileNameArray,1);
if FilePath
    if FileNum==1
%         set(handles.FileIndex,  'Enable', 'on');
        set(handles.FileIndex,  'Min', 0.5, 'Value', 1, 'Max', 1, 'SliderStep',[0,0]);
    elseif FileNum>1
%         set(handles.FileIndex,  'Enable', 'on');
        set(handles.FileIndex,  'Min', 1, 'Value', 1, 'Max', FileNum, 'SliderStep',[1/(FileNum-1),1/(FileNum-1)]);
    end
end
% get reference
Ref=0;
% if strcmpi(PrePar.SystemID, 'THQDOCT')
%     fid = fopen([FilePath,'RefAvrSpect.dat']);
%     Ref= fread(fid, 2048, 'double');
%     fclose(fid); 
% else
%     Ref = GetOCTData([FilePath,'fpnraw0000.raw'], PrePar.SystemID);
%     Ref = mean(Ref,2);
% end
function FileIndex_Callback(hObject, eventdata, handles)
global FilePath FileNameArray FileName PrePar IntPar;
global MIPImage;
% hObject    handle to FileIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileIndex=round(get(hObject,'Value'));
FileName=deblank(FileNameArray(FileIndex,:));
set(hObject,'Value', FileIndex);
set(handles.FileIndexText,'String',['FileIndex:',num2str(FileIndex)]);
set(handles.FileName,'String',FileName(1:end-4));
OCTData=GetOCTData([FilePath,FileName],PrePar.SystemID,1024, 1376); 
% OCTData=RemoveRef(OCTData,PrePar.ReferenceMode);
% split Spectrum SLO

% %%projection
% if mod(FileIndex, 2)
%     OCTData = rot90(OCTData,2);
% end
% M = 1;
% OCTData=RemoveRef(OCTData,'Mean');
% SubWidth = floor(size(OCTData,1)/M);
% if FileIndex == 1
%     MIPImage = zeros(200*M, 200);
% end
% for i = 1:M
%     MIPImage((0:199)*M+i,FileIndex) = sum(OCTData(i*SubWidth-SubWidth+1:i*SubWidth,:).^2);
% %     MIPImage(:,FileIndex,i) = MIPImage(:,FileIndex,i) /max(MIPImage(:,FileIndex,i));
% end
% %normalization
% if FileIndex == 200
%     NormMIPImage = MIPImage;
%     for i = 1:M
%         NormMIPImage(i:M:end,:) = NormMIPImage(i:M:end,:) / mean(mean(NormMIPImage(i:M:end,:)));
%     end   
%     NormMIPImage = imresize(NormMIPImage, [400, 200]);
%     % MIPImage(:,2:2:end) = flipud(MIPImage(:,2:2:end));
%     figure;imtool(circshift(NormMIPImage, [160, 0]));
% end 

if strcmpi(PrePar.SystemID, 'datfile')
    OCTData=imresize(OCTData,[1023 2048]);
    DisplayRange=GetDisplayRange(OCTData, IntPar.IntDispRange);
    figure;imshow(OCTData,'DisplayRange',DisplayRange);
    return;
else
    Reconstruct_Callback(handles, OCTData);
end
function Reconstruct_Callback(handles, RawData)
global PrePar IntPar DopplerPar ImgPar AgoPar CaliedData FilePath FileName ;
CaliedData=PreProcess(PrePar,  RawData, FilePath, FileName);
if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(CaliedData,PrePar,IntPar,ImgPar,FilePath,FileName,Const('Recon'));    
end
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(CaliedData,PrePar, DopplerPar,ImgPar,FilePath,FileName,Const('Recon'));
end
if get(handles.EnableAgo,'Value')
    ReconstructAngioFlow(CaliedData,PrePar, AgoPar,ImgPar,FilePath,FileName,Const('Recon'));
end
function BatchReconstruct_Callback(hObject, eventdata, handles)
FileNum=get(handles.FileIndex,'Max');
AutoReconstructProgressBar = waitbar(0,'Batch auto reconstruction...');
for i = get(handles.FileIndex,'Value'):get(handles.FileIndex,'Max')
    set(handles.FileIndex,  'Value',i);
    FileIndex_Callback(handles.FileIndex, eventdata, handles);
    waitbar(i/FileNum,AutoReconstructProgressBar);
end
close(AutoReconstructProgressBar);
function EnableIntensity_Callback(hObject, eventdata, handles)
function EnableDoppler_Callback(hObject, eventdata, handles)
function EnableAgo_Callback(hObject, eventdata, handles)
function GetParameters(handles)
global PrePar IntPar DopplerPar AgoPar ImgPar ;
PrePar.ReferenceMode=get(handles.ReferenceMode,'String');
PrePar.ReferenceMode=PrePar.ReferenceMode{get(handles.ReferenceMode,'Value')};
PrePar.GenDispCoeff=get(handles.GenDispCoeff,'Value');
PrePar.SystemID=get(handles.SystemID,'String');  
PrePar.SystemID=PrePar.SystemID{get(handles.SystemID,'Value')};
PrePar.DispFileName=get(handles.DispFileName,'String');
PrePar.DeConjugate=get(handles.DeConjugate,'String');
PrePar.DeConjugate=PrePar.DeConjugate{get(handles.DeConjugate,'Value')};
PrePar.DeConjOrder=get(handles.DeConjOrder,'String');
PrePar.DeConjOrder=str2num(PrePar.DeConjOrder(end));%% Image Processing
PrePar.DeConjWindow=get(handles.DeConjWindow,'String');
PrePar.DeConjWindow=str2num(PrePar.DeConjWindow(11:end));

IntPar.ImageScale=get(handles.IntScale,'String');
IntPar.ImageScale=IntPar.ImageScale{get(handles.IntScale,'Value')};
IntPar.IntDispRange=get(handles.IntDispRange,'String');
IntPar.IntDispRange=str2num(IntPar.IntDispRange(7:end));
IntPar.AutoSaveImg=get(handles.AutoSaveIntensityImage,'Value');
IntPar.AutoSaveFlt=get(handles.AutoSaveFloatData,'Value');


DopplerPar.Method=get(handles.DopplerMethod,'String');
DopplerPar.Method=DopplerPar.Method{get(handles.DopplerMethod,'Value')};
DopplerPar.IntThresh=str2double( get(handles.DopplerIntThresh,'String'));
DopplerPar.FilterType=get(handles.DopplerFilterType,'String');
DopplerPar.FilterType=DopplerPar.FilterType{get(handles.DopplerFilterType,'Value')};
DopplerPar.FilterSize=get(handles.DopplerFilterSize,'String');
DopplerPar.FilterSize=str2num(DopplerPar.FilterSize(12:end));
DopplerPar.RemoveBMA=get(handles.RemoveBMA,'Value');
DopplerPar.AbsFlow=get(handles.AbsFlow,'Value');
DopplerPar.AutoSaveImg=get(handles.AutoSaveDopplerFlow,'Value');
DopplerPar.AutoSaveFlt=get(handles.AutoSaveDopplerData,'Value');

AgoPar.Method=get(handles.AgoMethod,'String');
AgoPar.Method=AgoPar.Method{get(handles.AgoMethod,'Value')};
AgoPar.AutoSaveImg=get(handles.AutoSaveAgoImg,'Value');
AgoPar.AutoSaveFlt=get(handles.AutoSaveAgoData,'Value');

ImgPar.ImgFilterType=get(handles.ImgFilterType,'String');
ImgPar.ImgFilterType=ImgPar.ImgFilterType{get(handles.ImgFilterType,'Value')};
ImgPar.ImgFilterSize=get(handles.ImgFilterSize,'String');
ImgPar.ImgFilterSize=str2num(ImgPar.ImgFilterSize(12:end));
%ImgPar.fftShift=get(handles.fftShift,'Value');
ImgPar.ZeroPadNum=str2double(get(handles.ZeroPadNum,'String'));
ImgPar.ImageFOV=get(handles.ImageFOV,'String');
ImgPar.ImageFOV=str2num(ImgPar.ImageFOV(5:end));
ImgPar.ImgResize=get(handles.ImgResize,'String');
ImgPar.ImgResize=str2num(ImgPar.ImgResize(8:end));
%% PreProcess
function ReferenceMode_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.ReferenceMode=get(handles.ReferenceMode,'String');
PrePar.ReferenceMode=PrePar.ReferenceMode{get(handles.ReferenceMode,'Value')};
function GenDispCoeff_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.GenDispCoeff=get(handles.GenDispCoeff,'Value');
function DispFileName_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.DispFileName=get(handles.DispFileName,'String');
function DeConjugate_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.DeConjugate=get(handles.DeConjugate,'String');
PrePar.DeConjugate=PrePar.DeConjugate{get(handles.DeConjugate,'Value')};
function DeConjOrder_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.DeConjOrder=get(handles.DeConjOrder,'String');
PrePar.DeConjOrder=str2num(PrePar.DeConjOrder(end));%% Image Processing
function DeConjWindow_Callback(hObject, eventdata, handles)
global PrePar;
PrePar.DeConjWindow=get(handles.DeConjWindow,'String');
PrePar.DeConjWindow=str2num(PrePar.DeConjWindow(11:end));%% Image Processing
%% IntensityImage 
function IntScale_Callback(hObject, eventdata, handles)
global CaliedData  PrePar IntPar ImgPar FilePath FileName;
IntPar.ImageScale=get(handles.IntScale,'String');
IntPar.ImageScale=IntPar.ImageScale{get(handles.IntScale,'Value')};

if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(CaliedData,PrePar,IntPar,ImgPar,FilePath,FileName,Const('Recon'));    
end
function IntDispRange_Callback(hObject, eventdata, handles)
global IntPar  PrePar ImgPar FilePath FileName;
IntPar.IntDispRange=get(handles.IntDispRange,'String');
IntPar.IntDispRange=str2num(IntPar.IntDispRange(7:end));
if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(Const('Null'),PrePar, IntPar,ImgPar,FilePath,FileName,Const('Display'));    
end
function AutoSaveIntensityImage_Callback(hObject, eventdata, handles)
global IntPar;
IntPar.AutoSaveImg=get(handles.AutoSaveIntensityImage,'Value');
function AutoSaveFloatData_Callback(hObject, eventdata, handles)
global IntPar;
IntPar.AutoSaveFlt=get(handles.AutoSaveFloatData,'Value');
%% Doppler Flow functions
function DopplerMethod_Callback(hObject, eventdata, handles)
global CaliedData  PrePar DopplerPar ImgPar FilePath FileName;
DopplerPar.Method=get(handles.DopplerMethod,'String');
DopplerPar.Method=DopplerPar.Method{get(handles.DopplerMethod,'Value')};

if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(CaliedData,PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Recon'));
end
function DopplerFilterType_Callback(hObject, eventdata, handles)
global CaliedData  PrePar DopplerPar ImgPar FilePath FileName;
DopplerPar.FilterType=get(handles.DopplerFilterType,'String');
DopplerPar.FilterType=DopplerPar.FilterType{get(handles.DopplerFilterType,'Value')};
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(CaliedData,PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Recon'));
end
function DopplerFilterSize_Callback(hObject, eventdata, handles)
global CaliedData PrePar DopplerPar ImgPar FilePath FileName;
DopplerPar.FilterSize=get(handles.DopplerFilterSize,'String');
DopplerPar.FilterSize=str2num(DopplerPar.FilterSize(12:end));
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(CaliedData,PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Recon'));
end
function DopplerIntThresh_Callback(hObject, eventdata, handles)
global DopplerPar  PrePar ImgPar FilePath FileName;
DopplerPar.IntThresh=str2double( get(handles.DopplerIntThresh,'String'));
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(Const('Null'),PrePar,DopplerPar,ImgPar,FilePath, FileName,Const('Thresh'));
end
function RemoveBMA_Callback(hObject, eventdata, handles)
global CaliedData PrePar  DopplerPar ImgPar FileName;
DopplerPar.RemoveBMA=get(handles.RemoveBMA,'Value');

if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(CaliedData,PrePar,DopplerPar,ImgPar,FileName,Const('Recon'));
end
function AbsFlow_Callback(hObject, eventdata, handles)
global DopplerPar  PrePar ImgPar FileName;

DopplerPar.AbsFlow=get(handles.AbsFlow,'Value');
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(Const('Null'),PrePar,DopplerPar,ImgPar,FileName,Const('Filter'));
end
function AutoSaveDopplerFlow_Callback(hObject, eventdata, handles)
global DopplerPar;
DopplerPar.AutoSaveImg=get(handles.AutoSaveDopplerFlow,'Value');
function AutoSaveDopplerData_Callback(hObject, eventdata, handles)
global DopplerPar;
DopplerPar.AutoSaveFlt=get(handles.AutoSaveDopplerData,'Value');
%% Angio Functions
function AutoSaveAgoImg_Callback(hObject, eventdata, handles)
global AgoPar;
AgoPar.AutoSaveImg=get(handles.AutoSaveAgoImg,'Value');
function AutoSaveAgoData_Callback(hObject, eventdata, handles)
global AgoPar;
AgoPar.AutoSaveFlt=get(handles.AutoSaveAgoData,'Value');
function AgoMethod_Callback(hObject, eventdata, handles)
global AgoPar;
AgoPar.Method=get(handles.AgoMethod,'String');
AgoPar.Method=AgoPar.Method{get(handles.AgoMethod,'Value')};

%% Image enhance
function ImgFilterType_Callback(hObject, eventdata, handles)
global DopplerPar  PrePar IntPar ImgPar FilePath FileName;
ImgPar.ImgFilterType=get(handles.ImgFilterType,'String');
ImgPar.ImgFilterType=ImgPar.ImgFilterType{get(handles.ImgFilterType,'Value')};

if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(Const('Null'),PrePar, IntPar,ImgPar,FilePath, FileName,Const('Filter'));    
end
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(Const('Null'),PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Filter'));
end
function ImgFilterSize_Callback(hObject, eventdata, handles)
global DopplerPar  PrePar IntPar ImgPar FilePath FileName;
ImgPar.ImgFilterSize=get(handles.ImgFilterSize,'String');
ImgPar.ImgFilterSize=str2num(ImgPar.ImgFilterSize(12:end));
if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(Const('Null'),PrePar, IntPar,ImgPar,FilePath,FileName,Const('Filter'));    
end
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(Const('Null'),PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Filter'));
end
function ImgResize_Callback(hObject, eventdata, handles)
global ImgPar  PrePar IntPar DopplerPar FilePath FileName;
ImgPar.ImgResize=get(handles.ImgResize,'String');
ImgPar.ImgResize=str2num(ImgPar.ImgResize(8:end));
if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(Const('Null'),PrePar, IntPar,ImgPar,FilePath,FileName,Const('Filter'));    
end
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(Const('Null'),PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Filter'));
end
function ImageFOV_Callback(hObject, eventdata, handles)
global ImgPar  PrePar IntPar DopplerPar FilePath FileName;
ImgPar.ImageFOV=get(handles.ImageFOV,'String');
ImgPar.ImageFOV=str2num(ImgPar.ImageFOV(5:end));
if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(Const('Null'),PrePar, IntPar,ImgPar,FilePath,FileName,Const('Filter'));    
end
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(Const('Null'),PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Thresh'));
end
function ZeroPadNum_Callback(hObject, eventdata, handles)
global CaliedData  PrePar ImgPar IntPar DopplerPar FilePath  FileName;
ImgPar.ZeroPadNum=str2double(get(handles.ZeroPadNum,'String'));
if get(handles.EnableIntensity,'Value')
    ReconstructIntensityImage(CaliedData,PrePar,IntPar,ImgPar,FilePath,FileName,Const('Recon'));    
end
if get(handles.EnableDoppler,'Value')
    ReconstructDopplerFlow(CaliedData,PrePar,DopplerPar,ImgPar,FilePath,FileName,Const('Recon'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Misc Tools                                    % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CombineBFFrames_Callback(hObject, eventdata, handles)
global PrePar;
CombineBackForwardFrames(PrePar.SystemID);
function LoadCalibration_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCalibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CaliCoeff;
[FileName,FilePath]=uigetfile({'*.mat';'*.*'}, 'Pick calibration data file');
if FilePath
    if strcmpi(FileName(end-2:end),'mat')
        load([FilePath,FileName]);
    else
        GenCaliCoeff(load([FilePath,FileName]));
    end
else
    CaliCoeff=[];
end

