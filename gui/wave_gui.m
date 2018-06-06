function varargout = wave_gui(varargin)
% WAVE_GUI MATLAB code for wave_gui.fig
%      WAVE_GUI, by itself, creates a new WAVE_GUI or raises the existing
%      singleton*.
%
%      H = WAVE_GUI returns the handle to a new WAVE_GUI or the handle to
%      the existing singleton*.
%
%      WAVE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVE_GUI.M with the given input arguments.
%
%      WAVE_GUI('Property','Value',...) creates a new WAVE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wave_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wave_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wave_gui

% Last Modified by GUIDE v2.5 01-Jun-2018 13:41:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wave_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @wave_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before wave_gui is made visible.
function wave_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wave_gui (see VARARGIN)

% Choose default command line output for wave_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wave_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function [wave_data]=pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(lower(handles.edit1.String),'wave csv file')
    [file, path] = uigetfile({'*.csv';'*.*'},'File Selector');
    % edit1_Callback(hObject, eventdata, handles)
    handles.edit1.String=sprintf('%s/%s',path,file);
end

% [fname, path] = uigetfile({'*.csv';'*.*'},'File Selector');
% filename = set(hObject,'String',sprintf('%s/%s',path,fname));
filename = handles.edit1.String;

% filename = sprintf('%s/tek%04.0fCH1.csv',fpath,n);
% [~,fname,fext] = fileparts(filename);
% outfname = sprintf('%s/%s.xlsx',fpath,fname);

fid=fopen(sprintf('%s',filename),'r');
if fid<1; fprintf('%s could not be found',filename); return; end

%initialize variables used in readin the file

header_flag = true;
header = '';

% read the file
tline = fgetl(fid);
while ischar(tline)    
    tline = fgetl(fid);
   
    %Check for end of header header
    if strcmp(tline,'TIME,CH1'); header_flag = false; end
    
    % read the data appropriately based on header or not
    if header_flag
        header = sprintf('%s\n%s',header,tline);
    else
        time_voltage_data=textscan(fid,'%f,%f');
        time = time_voltage_data{1};
        ch1 = time_voltage_data{2};
        break
    end
end


% time
time(isinf(ch1))=[]; %remove unreal values
time_micro = time*1e6; % Time in microseconds

% voltage
ch1(isinf(ch1))=[]; %remove unreal values
voltage_mean = mean(ch1);
voltage_relative = ch1 - voltage_mean;

% sensitivity
calibration = str2double(handles.edit2.String); %(V/MPa)

% Pressure = ch1/calibration*1e6;
Pressure_MPa = voltage_relative/calibration;
Pressure = Pressure_MPa*1e6;


% tabplot('Pressure')
axes(handles.axes1)

plot(time_micro,Pressure_MPa)
xlim([min(time_micro), max(time_micro)])
xlabel('time (\mus)')
ylabel('Pressure (MPa)')

rho = 1e3;
c=1.5e3;
Velocity = Pressure / (rho*c);
Instant_Intensity = Pressure.*Velocity;
Instant_Intensity_Wcm2 = Instant_Intensity * (0.01)^2;
Intensity = Pressure.^2/(rho*c);
PII = trapz(time, Instant_Intensity)/range(time); %Pulse intensity integral based on Kinsler Eq. (5.9.1)
PII_cm = trapz(time, Instant_Intensity)/range(time)/1e4; %Pulse intensity integral based on Kinsler Eq. (5.9.1) per centimeter


date_rec = '02-May-2018';

tekfile = filename;

sample_duration_micro = range(time_micro);
% pulse_duration_micro =

Pmax_MPa = max(Pressure_MPa);
Pmin_MPa = min(Pressure_MPa);
Pabs_mid_MPa = (Pmax_MPa - Pmin_MPa)/2;

wavedata = struct('time',time,'time_micro',time*1e6,'voltage','voltage_mean',voltage_mean,ch1,'calibration_VoltMPa',calibration,...
    'Pressure',Pressure,'Pressure_MPa',Pressure_MPa,'PressurePeakPositive',max(Pressure),'PressurePeakNegative',min(Pressure),'PressureCenter',...
    'Intensity',Instant_Intensity,'Intensity_Wcm2',Instant_Intensity_Wcm2,'IntensityPeak_Wcm2',max(Instant_Intensity_Wcm2), 'IntensityPulseAverage_Wcm2',0,...
    'WaveDuration',range(time),'WaveDuration_micro',range(time_micro),'PulseDuration',0,'PulseDuration_micro',0,...
    'Pressure_omega',[],'frequency_spectrum',[],...
    'tekfile',filename...
    );

hObject.UserData = wavedata;


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=1+1;


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,~]=ginput(1);
x = x(end);
xstr=sprintf('%.3f',x);
handles.edit4.String = xstr;
axes(handles.axes1); hold on;
if isfield(handles.axes1.UserData,'LeftBoundaryPlot'); delete(handles.axes1.UserData.LeftBoundaryPlot); end % remove old line
handles.axes1.UserData.LeftBoundaryPlot = plot(x*[1,1], handles.axes1.YLim,'r--'); hold off;
a=1+1;


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB




% handles    structure with handles and user data (see GUIDATA)
[x,~]=ginput(1);
x = x(end);
xstr=sprintf('%.3f',x);
handles.edit3.String = xstr;
axes(handles.axes1); hold on;
plot(x*[1,1], handles.axes1.YLim,'r--'); hold off


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(lower(handles.edit1.String),'wave csv file')
    [file, path] = uigetfile({'*.csv';'*.*'},'File Selector');
    % edit1_Callback(hObject, eventdata, handles)
    handles.edit1.String=sprintf('%s/%s',path,file);
end


















function [ freqdata ] = bfft( timedata )
    [rows,columns]=size(timedata);
    
    if columns>rows %I want each column to be the complete fourier information
        timedata=transpose(timedata);
    end
    
    L=round(max(size(timedata))/2)*2; %Account for possibility that L is not a power of 2
    NN=min(size(timedata));
    freqdata=zeros(L/2+1,NN);
    
    for i=1:NN
        temptimedata=timedata(:,i); %Assign the ith column
        temp=conj(fft(temptimedata)); %FFT Convention enforced
        freqdata(:,i)=temp(1:L/2+1); %Grab only the information before the folding frequency
    end
    
    if columns>rows %switch it back if necessary
        freqdata=transpose(freqdata);
    end





function [ timedata ] = bifft( freqdata )
    [rows,columns]=size(freqdata);
    
    if columns>rows %I want each column to be the complete fourier information
        freqdata=transpose(freqdata);
    end
    
    L=2*max(size(freqdata))-2;
    NN=min(size(freqdata));
    timedata=zeros(L,NN);
    
    for i=1:NN
        tempfreqdata=freqdata(:,i); %Assign the ith column
        
        tempfreqdata(1)=real(tempfreqdata(1)); %Force the DC bin to be real, otherwise the time domain has a constant imaginary component. As below, a sign change here doesn't change the realness. Not sure what's "most" correct lol...

        L=2*length(tempfreqdata)-2; %Length of time output, such that Length(freq)=L/2+1
        preIFFT=zeros(L,1); %Prepare a full-length vector for IFFT
        preIFFT(1:L/2+1)=tempfreqdata; %Assign the freq data to first half
        temp=conj(flipud(tempfreqdata)); %Flip freq data, conjugate and store
        preIFFT(L/2+1:L)=temp(1:L/2); %Assign the flipped data, skipping the first data point
        preIFFT(L/2+1)=real(preIFFT(L/2+1)); %Force the folding frequency to be real (abs or real part?) - I'm choosing the real part, since having a negative number here doesn't destroy the real-ness of the time data
        timedata(:,i)=ifft(conj(preIFFT)); %Fixed for FFT Convention
    end
    
    if columns>rows %switch it back if necessary
        timedata=transpose(timedata);
    end

