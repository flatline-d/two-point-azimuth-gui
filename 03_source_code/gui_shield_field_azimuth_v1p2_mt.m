function varargout = gui_shield_field_azimuth_v1p2_mt(varargin)
% GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT MATLAB code for gui_shield_field_azimuth_v1p2_mt.fig
%      GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT, by itself, creates a new GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT or raises the existing
%      singleton*.
%
%      H = GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT returns the handle to a new GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT or the handle to
%      the existing singleton*.
%
%      GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT.M with the given input arguments.
%
%      GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT('Property','Value',...) creates a new GUI_SHIELD_FIELD_AZIMUTH_V1P2_MT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_shield_field_azimuth_v1p2_mt_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_shield_field_azimuth_v1p2_mt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_shield_field_azimuth_v1p2_mt

% Last Modified by GUIDE v2.5 14-Apr-2016 15:35:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_shield_field_azimuth_v1p2_mt_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_shield_field_azimuth_v1p2_mt_OutputFcn, ...
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


% --- Executes just before gui_shield_field_azimuth_v1p2_mt is made visible.
function gui_shield_field_azimuth_v1p2_mt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_shield_field_azimuth_v1p2_mt (see VARARGIN)

% Choose default command line output for gui_shield_field_azimuth_v1p2_mt
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_shield_field_azimuth_v1p2_mt wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_shield_field_azimuth_v1p2_mt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_data_selection.
function pushbutton_data_selection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_data_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Use "uigetfile" to activate user selection window.
% [FileName,PathName} = uigetfile displays a modal dialog box lists files
% in the current folder and enables the user to select a file (or enter its
% name).
[fileName_shield_field_pts, pathname1] = uigetfile({'*.csv';'*.txt'},'Select shield field location file')
%
% Assign this value to display in the string field of this text string
set(handles.editField_shield_field_file_name, 'string', fullfile(pathname1, fileName_shield_field_pts));
%
% Execute (or reexecute) callback to text field display to update with
% user-selected file name.
editField_shield_field_file_name_Callback(hObject, eventdata, handles)


function editField_shield_field_file_name_Callback(hObject, eventdata, handles)
% hObject    handle to editField_shield_field_file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editField_shield_field_file_name as text
%        str2double(get(hObject,'String')) returns contents of editField_shield_field_file_name as a double


% --- Executes during object creation, after setting all properties.
function editField_shield_field_file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editField_shield_field_file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_shields.
function pushbutton_plot_shields_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_shields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
axes(handles.axes1_xy_locations); % create a graphical object
cla; % clear current axes
%
% Retrieve previously selected shield field position data file
sf_xy_pts_name=get(handles.editField_shield_field_file_name, 'String');
%
% Call function "plot_xy_point_data.m" to read in the data file, plot
% the shield locations in a scatterplot, and return the number of shields
npoints = plot_xy_point_data(sf_xy_pts_name)
%
% Assign # shields to display in the string field of this text string
set(handles.plot_label1,'String', ['Number of shields = ', num2str(npoints)])


% --- Executes on button press in pushbutton_compute_az.
function pushbutton_compute_az_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_compute_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Update status string to read 'busy' during execution
set(handles.status_text,'String','Status: busy');
pause(0.1)
%
% Retrieve previously selected shield field position data file
sf_xy_pts_name=get(handles.editField_shield_field_file_name, 'String');
% Retrieve planetary body selection
planet=get(handles.text_planet_selection, 'String');
% Retrieve azimuth method selection (Lutz or Cebria method)
% Retrieves the Tag property of the selected radio button in buttongroup1
% and assigns it to variable "method" as a string.
method=get(get(handles.uibuttongroup1,'SelectedObject'), 'Tag')
%
% Use switch statement to determine which radiobutton was selected and
% execute appropriate azimuth method.
switch method
    case 'radiobutton6_lutz'
        % Designate which set of axes are active (i.e., will have plots added to
        % them)
        azh1 = handles.axes3;
        azh2 = handles.axes4;
        % Call function "azimuth_method_cebria_v3_mt.m" to tabulate the azimuth values
        % and create two plots: a histogram in azh1 and a rose plot in azh2.
        % Requires Mapping Toolbox.
        [az_h_vals, bin_centers] = azimuth_method_v3_mt(sf_xy_pts_name, azh1, azh2, planet);
    case 'radiobutton7_cebria'
        azh1 = handles.axes1_xy_locations;
        azh2 = handles.axes3;
        azh3 = handles.axes4;
        % Call function "azimuth_method_cebria.m" to tabulate the azimuth values
        % and create three plots: (1) line segments in x-y location plot in azh1,
        % (2) a histogram in azh2, and (3) a rose plot in azh3.
        [az_h_vals, bin_centers] = azimuth_method_cebria_mt(sf_xy_pts_name, azh1, azh2, azh3, planet);
        %
        % Overlay x-y location plot over line segment plot in azh1.
        axes(handles.axes1_xy_locations); % create a graphical object
        hold on;
        % Call function "plot_xy_point_data.m" to read in the data file, plot
        % the shield locations in a scatterplot, and return the number of shields
        npoints = plot_xy_point_data(sf_xy_pts_name)
end
%
% Store observed azimuth histogram "az_h_vals" in the application workspace
% using the SETAPPDATA function. "figure1" is tag of GUI object
setappdata(handles.figure1, 'observed_hist', az_h_vals)
setappdata(handles.figure1, 'observed_hist_bins', bin_centers)
%
% Reset status string value
pause(0.1)
set(handles.status_text,'string','Status: done/paused');


% --- Executes during object creation, after setting all properties.
function pushbutton_compute_az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_compute_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4


% --- Executes during object creation, after setting all properties.
function plot_label1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_label1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function status_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_start_monte_carlo.
function pushbutton_start_monte_carlo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start_monte_carlo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Update status string to read 'busy' during execution
set(handles.status_text,'String','Status: busy');
pause(0.1)
%
% Retrieve previously selected shield field position data file
sf_xy_pts_name=get(handles.editField_shield_field_file_name, 'String');
%
% Retrieve planetary body selection
planet=get(handles.text_planet_selection, 'String');
%
% Run azimuth pushbutton as a failsafe / bug prevention
pushbutton_compute_az_Callback(handles.pushbutton_compute_az, eventdata, handles)
%
% Retrieve previously calcualted histogram of observed data
obs_az_hist = getappdata(handles.figure1, 'observed_hist');
%
% Retrieve number of Monte Carlo runs to perform from exitbox associated
% with slider, convert from string to a number
num_MC_runs = str2num(get(handles.editbox_number_MC, 'String'));
%
% Designate which set of axes are active (i.e., will have plots added to
% them)
azh7 = handles.axes7;
%
% Retrieve azimuth method selection (Lutz or Cebria method)
% Retrieves the Tag property of the selected radio button in buttongroup1
% and assigns it to variable "method" as a string.
method=get(get(handles.uibuttongroup1,'SelectedObject'), 'Tag')
%
% Use switch statement to determine which radiobutton was selected and
% execute appropriate azimuth method.
switch method
    case 'radiobutton6_lutz'
        % Call function "monte_carlo_model_v4_mt.m" to run through a series of Monte
        % Carlo models (requires Mapping Toolbox).
        [bn,nh,th] = monte_carlo_model_v4_mt(sf_xy_pts_name, obs_az_hist, azh7, num_MC_runs, planet);
    case 'radiobutton7_cebria'
        % Call function "monte_carlo_model_v3_cebria_mt" to run through a series of Monte
        % Carlo models.
        [bn,nh,th] = monte_carlo_model_v4_cebria_mt(sf_xy_pts_name, obs_az_hist, azh7, num_MC_runs, planet);
end
%
% Store Monte Carlo model output in the application workspace
% using the SETAPPDATA function. "figure1" is tag of GUI object
setappdata(handles.figure1, 'MC_bins', bn); % histogram bin center values
setappdata(handles.figure1, 'MC_norm_hist', nh); % normalized histogram
setappdata(handles.figure1, 'MC_thresh', th); % 95% threshold values
%
% Reset status string value
pause(0.1)
set(handles.status_text,'string','Status: done/paused');


% --- Executes on slider movement.
function slider_no_MC_runs_Callback(hObject, eventdata, handles)
% hObject    handle to slider_no_MC_runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%
% Get value of slider. Note slider is exponential between 1 and 1000 (i.e.,
% 10^0 and 10^3). Value of the slider is taken as base 10 exponent to
% assign value to the edit text field.
MC_runs_exponent = get(handles.slider_no_MC_runs, 'Value');
MC_runs = 10^(MC_runs_exponent);
%
% Set retrieved value to edit text field adjacent to slider
set(handles.editbox_number_MC, 'String', sprintf('%.0f',MC_runs));
%
% Execute (or reexecute) callback to text field display to update with
% user-selected file name.
editbox_number_MC_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider_no_MC_runs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_no_MC_runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editbox_number_MC_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_number_MC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbox_number_MC as text
%        str2double(get(hObject,'String')) returns contents of editbox_number_MC as a double


% --- Executes during object creation, after setting all properties.
function editbox_number_MC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editbox_number_MC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_saveMC.
function pushbutton_saveMC_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveMC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% This push button save the Monte Carlo model result to a comma-delimited
% text file and detached header with column titles
%
% Retreive data stored in the application workspace
% using the GETAPPDATA function. "figure1" is tag of GUI object
bn = getappdata(handles.figure1, 'MC_bins'); % histogram bin center values
nh = getappdata(handles.figure1, 'MC_norm_hist'); % normalized histogram
th = getappdata(handles.figure1, 'MC_thresh'); % 95% threshold values
M2 = cat(2, bn', nh, th); % concatenate 3 arrays into matrix M2
%
% Query the user for the file name and location of the to-be-saved file 
[FileName2,PathName2] = uiputfile('*.csv','Save Monte Carlo results','MC_hist');
[p,n2,e] = fileparts(FileName2); % split file name into base name and
                               % extension; p = path, n=name, e=extension
                               % path 'p' is empty
%
% Output data as a text file with deteached header
% 'dlmwrite' writes matrix data to delimited ascii file, default delimiter
% is a comma.
FileFullPath2 = [PathName2 FileName2];  % full_path + file_name.csv
dlmwrite(FileFullPath2, M2);  
%
% Write column headers of ascii data file to separate label file
% HdrSaveName2 = full_path + file_name + _header.txt
HdrSaveName2 = [PathName2 n2 '_header.txt']; 
fid = fopen(HdrSaveName2, 'wt');
fprintf(fid, 'Bin_center_(degrees), Normalized_Num_azimuths, threshold_value \n');
fclose(fid)


% --- Executes on button press in pushbutton_saveHist.
function pushbutton_saveHist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% This push button save the raw azimuth histogram to a text file
%
% Retreive azimuth histogram data stored in the application workspace
% using the GETAPPDATA function. "figure1" is tag of GUI object
az_hist = getappdata(handles.figure1, 'observed_hist');
bin_centers = getappdata(handles.figure1, 'observed_hist_bins');
M = cat(2, bin_centers', az_hist); % concatenate 2 arrays into matrix M
%
% Query the user for the file name and location of the to-be-saved file 
[FileName,PathName] = uiputfile('*.csv','Save azimuth histogram','az_hist');
[p,n,e] = fileparts(FileName); % split file name into base name and
                               % extension; p = path, n=name, e=extension
                               % path 'p' is empty
%
% Output data as a text file with deteached header
% 'dlmwrite' writes matrix data to delimited ascii file, default delimiter
% is a comma.
FileFullPath = [PathName FileName]; % full_path + file_name.csv
dlmwrite(FileFullPath, M);  
%
% Write column headers of ascii data file to separate label file
HdrSaveName = [PathName n '_header.txt']; % full_path + file_name + _header.txt
fid = fopen(HdrSaveName, 'wt');
fprintf(fid, 'Bin_center_(degrees), Number_of_azimuths \n');
fclose(fid)



% --- Executes during object creation, after setting all properties.
function radiobutton1_venus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1_venus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton2_earth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2_earth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton3_mars_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton3_mars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton6_lutz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton6_lutz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton7_cebria_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton7_cebria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel6_choose_planet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6_choose_planet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel6_choose_planet.
function uipanel6_choose_planet_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6_choose_planet 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%
if (hObject == handles.radiobutton8_unitsphere)
    set(handles.text_planet_selection, 'String', 'Unit Sphere');
    %
elseif (hObject == handles.radiobutton1_venus)
    set(handles.text_planet_selection, 'String', 'Venus');
    %
elseif (hObject == handles.radiobutton2_earth)
    set(handles.text_planet_selection, 'String', 'Earth');
    %
elseif (hObject == handles.radiobutton3_mars)
    set(handles.text_planet_selection, 'String', 'Mars');
    %
else
    set(handles.text_planet_selection, 'String', 'Bad selection');
end


function text_planet_selection_Callback(hObject, eventdata, handles)
% hObject    handle to text_planet_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_planet_selection as text
%        str2double(get(hObject,'String')) returns contents of text_planet_selection as a double


% --- Executes during object creation, after setting all properties.
function text_planet_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_planet_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1_venus.
function radiobutton1_venus_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1_venus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1_venus


% --- Executes on button press in radiobutton6_lutz.
function radiobutton6_lutz_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6_lutz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of radiobutton6_lutz
%
% Get all the handles to everything we want to re-enable set into a single
% array, undoing disablement that happens if Cebria radio button is
% selected.
handlesArray = [handles.pushbutton_start_monte_carlo, handles.pushbutton_saveMC];
% Set the Enable parameter to on (i.e., re-enable them)
set(handlesArray, 'Enable', 'on');
% Also make the Monte Carlo plot visibl (axes7)
set(handles.axes7,'visible','on'), 


% --- Executes on button press in radiobutton7_cebria.
function radiobutton7_cebria_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7_cebria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of radiobutton7_cebria
%
% % Due to issues with applying the Monte Carlo results to the Cebria method,
% % disable the Monte Carlo pushbuttons when the Cebria radio button is
% % selected.
% % Get all the handles to everything we want to disable set into a single array.
% handlesArray = [handles.pushbutton_start_monte_carlo, handles.pushbutton_saveMC];
% % Set them all disabled (i.e., not enabled)
% set(handlesArray, 'Enable', 'off');
% % Also hide the Monte Carlo plot (axes7)
% set(handles.axes7,'visible','off'), 



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton1_venus.
function radiobutton1_venus_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1_venus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton6_lutz'
        display('Lutz method');
    case 'radiobutton7_cebria'
        display('Cebria method');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_planet_selection.
function text_planet_selection_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text_planet_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uibuttongroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton_data_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_data_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in radiobutton8_unitsphere.
function radiobutton8_unitsphere_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8_unitsphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8_unitsphere


% --- Executes during object creation, after setting all properties.
function axes7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes7
