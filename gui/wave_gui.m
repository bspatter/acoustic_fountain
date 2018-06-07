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

% Last Modified by GUIDE v2.5 06-Jun-2018 15:08:27

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

% clean-up anything we might not want
clearstuff(handles)

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



axes(handles.axes1)
handles.axes1.UserData.MainPlot = plot(time_micro,Pressure_MPa);
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

% FFT stuff
time_adjusted = time-time(1); % start time from 0;
TimeTotal=range(time_adjusted); % Total time
SamplingFrequency = 1/mean(diff(time_adjusted)); % sampling frequency
N = round(TimeTotal*SamplingFrequency); %number of data points
N = round(N/2)*2; %force it to be even
TimeTotal = N/SamplingFrequency; %new total time domain length, if N changed
TimeFFT = (0:(N-1))/SamplingFrequency; % time vector (max(t) is 1 bin short of T)
FrequencySpectrum = (0:(N/2))/N*SamplingFrequency; %frequency vector (min(f)=0, the DC component. max(f) = fs/2 exactly)
% Match pressure to output variables
if length(TimeFFT)==length(Pressure)
    y = Pressure;%function(t); %hypothetical time domain vector, length(y)=N
else
    y = interp1(time,Pressure, TimeFFT);%function(t); %hypothetical time domain vector, length(y)=N
end
PressureFFT= bfft(y);
[~,dominant_frequency_index]=max(abs(PressureFFT));
dominant_frequency = FrequencySpectrum(dominant_frequency_index);

% Is the signal undersampled?
undersampled_flag = true;
MI = 0;
if ~(max(FrequencySpectrum)==dominant_frequency)
    undersampled_flag = false;
%     MI = abs(min(pulse_pressure_MPa))/sqrt(pulse_frequency/1e6);
end

[Pressure_derated, time_derated] = derate_wave(time,Pressure,str2double(handles.edit7.String), str2double(handles.edit8.String));
Intensity_derated = Pressure_derated.^2/(rho*c);

wavedata = struct('time',time,'time_micro',time*1e6,'voltage',ch1,'voltage_mean',voltage_mean,'calibration_VoltMPa',calibration,...
    'Pressure',Pressure,'Pressure_MPa',Pressure_MPa,'PressurePeakPositive',max(Pressure),'PressurePeakNegative',min(Pressure),'PressureCenter',Pabs_mid_MPa,...
    'Intensity',Instant_Intensity,'Intensity_Wcm2',Instant_Intensity_Wcm2,'IntensityPeak_Wcm2',max(Instant_Intensity_Wcm2), 'IntensityPulseAverage_Wcm2',0,...
    'WaveDuration',range(time),'WaveDuration_micro',range(time_micro),'PulseDuration',0,'PulseDuration_micro',0,...
    'PressureFFT',PressureFFT(:),'FrequencySpectrum',FrequencySpectrum(:),'UnderSampled',undersampled_flag,...
    'Pressure_derated',Pressure_derated,'time_derated',time_derated,'Intensity_derated',Intensity_derated,...
    'rho',rho,'c',c,...
    'MI',MI,'tekfile',filename...
    );

handles.pushbutton1.UserData = wavedata;


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
wavedata = handles.pushbutton1.UserData;

pulse_interval = [str2double(handles.edit4.String), str2double(handles.edit3.String)]/1e6;
pulse_left_boundary_index = find(wavedata.time < pulse_interval(1),1,'last');
pulse_right_boundary_index = find(wavedata.time > pulse_interval(2),1,'first');

pulse_range_indeces = [pulse_left_boundary_index, pulse_right_boundary_index];
        
pulse_duration = range(wavedata.time(pulse_range_indeces));
pulse_duration_micro = pulse_duration*1e6;

pulse_pressure = wavedata.Pressure(pulse_range_indeces(1):pulse_range_indeces(2));
pulse_pressure_MPa = pulse_pressure/1e6;
pulse_intensity = wavedata.Intensity(pulse_range_indeces(1):pulse_range_indeces(2));
pulse_time = wavedata.time(pulse_range_indeces(1):pulse_range_indeces(2));
Half_Max_Intensity = 0.5*max(pulse_intensity);
pulse_half_times = intersections(pulse_time,(pulse_intensity-Half_Max_Intensity),[min(pulse_time),max(pulse_time)],[0,0]);
pulse_half_times_micro = pulse_half_times * 1e6;
pulse_half_intensities = interp1q(pulse_time,pulse_intensity,pulse_half_times);
pulse_half_intensities_cm = pulse_half_intensities/1e4;
FWHM = range(pulse_half_times);

pulse_average_intensity = trapz(pulse_time,pulse_intensity)/range(pulse_time);
pulse_average_intensity_cm = pulse_average_intensity / 1e4;

% FFT stuff to get frequency for the pulse alone
YP = fft(pulse_pressure);
L=length(pulse_pressure);
YP2 = abs(YP/L);
YP1 = YP2(1:(floor(L/2)+1));
YP1(2:end-1) = 2*YP1(2:end-1);

sampling_frequency = 1/mean(diff(wavedata.time));
frequency_spectrum = sampling_frequency*(0:(L/2))/L;
pulse_frequency = frequency_spectrum(YP1 == max(YP1));
if ~(max(frequency_spectrum)==pulse_frequency)
    MI = abs(min(pulse_pressure_MPa))/sqrt(pulse_frequency/1e6);
    fpulse = pulse_frequency;
else
    MI = 0;
    fpulse = 0;
end

% Assign relevant data to table in gui
pulse_table_data = [pulse_duration_micro, round(max(pulse_pressure_MPa),2), round(min(pulse_pressure_MPa),2), round((max(pulse_pressure_MPa) - min(pulse_pressure_MPa))/2,2), round(pulse_average_intensity_cm,1),round(fpulse/1e6,2),round(MI,2)]';
handles.uitable2.Data = pulse_table_data;
% Set width and height
handles.uitable2.ColumnWidth={125};
handles.uitable2.Position(3) = handles.uitable2.Extent(3);
handles.uitable2.Position(4) = handles.uitable2.Extent(4);


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
[file, path] = uigetfile({'*.csv';'*.*'},'File Selector');
% edit1_Callback(hObject, eventdata, handles)
handles.edit1.String=sprintf('%s/%s',path,file);


% --- Executes during object creation, after setting all properties.
function uibuttongroup2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
% function radiobutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wavedata = handles.pushbutton1.UserData;
switch (get(eventdata.NewValue,'Tag'))
    case 'radiobutton2'
%         handles.axes1.UserData.MainPlot = plot(wavedata.time_micro,wavedata.Pressure_MPa);
        handles.axes1.UserData.MainPlot.XData = wavedata.time_micro;
        handles.axes1.UserData.MainPlot.YData = wavedata.Pressure_MPa;
        handles.axes1.XLabel = xlabel('Time (\mus)');
        handles.axes1.YLabel = ylabel('Pressure (MPa)');
        handles.axes1.XLim = [min(wavedata.time_micro), max(wavedata.time_micro)];
        
        if handles.checkbox1.Value            
%             hold on
%             derated_plot = plot(wavedata.time_derated*1e6,wavedata.Pressure_derated/1e6,'r'); hold off
            handles.axes1.UserData.DeratedPlot.XData = wavedata.time_derated*1e6;
            handles.axes1.UserData.DeratedPlot.YData = wavedata.Pressure_derated/1e6;
%             handles.axes1.UserData.DeratedPlot = derated_plot;
        end
        
        
        
    case 'radiobutton3'
        if isfield(handles.axes1.UserData,'mytext'); delete(handles.axes1.UserData.mytext); end % Clear out unnecessary text
        axes(handles.axes1);
        handles.axes1.UserData.MainPlot.XData = wavedata.time_micro;
        handles.axes1.UserData.MainPlot.YData = wavedata.Intensity_Wcm2;
        handles.axes1.XLabel = xlabel('Time (\mus)');
        handles.axes1.YLabel =ylabel('Intensity (W / cm^2)');
        handles.axes1.XLim = [min(wavedata.time_micro), max(wavedata.time_micro)];
        
        if handles.checkbox1.Value            
%             hold on
%             derated_plot = plot(wavedata.time_derated*1e6,wavedata.Intensity_derated/1e4,'r'); hold off
            handles.axes1.UserData.DeratedPlot.XData = wavedata.time_derated*1e6;
            handles.axes1.UserData.DeratedPlot.YData = wavedata.Intensity_derated/1e4;
%             handles.axes1.UserData.DeratedPlot = derated_plot;
        end
        
        
    case 'radiobutton4'
        if isfield(handles.axes1.UserData,'mytext'); delete(handles.axes1.UserData.mytext); end % Clear out unnecessary text
        axes(handles.axes1);
        handles.axes1.UserData.MainPlot = plot(wavedata.FrequencySpectrum,abs(wavedata.PressureFFT));
%         handles.axes1.UserData.MainPlot.XData = wavedata.FrequencySpectrum;
%         handles.axes1.UserData.MainPlot.YData = wavedata.PressureFFT;
        handles.axes1.XLabel = xlabel('Frequency (Hz)');
        handles.axes1.YLabel =ylabel('FFT(Pressure) (Pa)');
        handles.axes1.XLim = [min(wavedata.FrequencySpectrum), max(wavedata.FrequencySpectrum)];
        if wavedata.UnderSampled
           handles.axes1.UserData.mytext = text(0.05,0.95,'Warning: Undersampled Signal','Units','Normalized','Color','red','FontWeight','bold');
        end
        
end

handles.axes1.YLimMode = 'auto';

function clearstuff(handles)
    if isfield(handles.axes1.UserData,'mytext'); delete(handles.axes1.UserData.mytext); end % Clear out unnecessary text


% --- Executes during object creation, after setting all properties.
function uitable2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
wavedata = handles.pushbutton1.UserData;

[wavedata.Pressure_derated, wavedata.time_derated] = derate_wave(wavedata.time,wavedata.Pressure,str2double(handles.edit7.String), str2double(handles.edit8.String));
wavedata.Intensity_derated = wavedata.Pressure_derated.^2/(wavedata.rho*wavedata.c);
handles.pushbutton1.UserData = wavedata;


if eventdata.Source.Value 
    switch handles.uibuttongroup2.SelectedObject.Tag
        case 'radiobutton2'
            axes(handles.axes1); hold on
%             derated_plot = plot(wavedata.time_derated, wavedata.Pressure_derated,'r'); hold off
            handles.axes1.UserData.DeratedPlot=plot(wavedata.time_derated*1e6, wavedata.Pressure_derated/1e6,'r'); hold off;

        case 'radiobutton3'
            axes(handles.axes1); hold on
            derated_plot = plot(wavedata.time_derated*1e6,wavedata.Intensity_derated/1e4,'r'); hold off
            handles.axes1.UserData.DeratedPlot = derated_plot;

        otherwise
    end
else
    if isfield(handles.axes1.UserData, 'DeratedPlot')
        delete(handles.axes1.UserData.DeratedPlot);
    end
end
handles.axes1.YLimMode = 'auto';

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
