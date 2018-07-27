% *******************************************************************
% *                        WAVEMAKER                                *
% * Software to simulate falling films using the Integral boundary
% * layer method ....
% *
% *                                                                 *
% *                   Wilko Rohlfs, Manuel Rietz                    *
% *  Institute of Heat and Mass Transfer, RWTH Aachen University    *
% *                       and Benoit Scheid
% *       TIPs lab – Microfluidics, Université Libre de Bruxelles   *
% *******************************************************************
%
% Copyright 2017 Wilko Rohlfs, Manuel Rietz and Benoit Scheid
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% This script contains the main graphical user interface (GUI) of the
% software Wavemaker

function varargout = WaveMaker(varargin)
% WAVEMAKER MATLAB code for WaveMaker.fig
%      WAVEMAKER, by itself, creates a new WAVEMAKER or raises the existing
%      singleton*.
%
%      H = WAVEMAKER returns the handle to a new WAVEMAKER or the handle to
%      the existing singleton*.
%
%      WAVEMAKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVEMAKER.M with the given input arguments.
%
%      WAVEMAKER('Property','Value',...) creates a new WAVEMAKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WaveMaker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WaveMaker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WaveMaker

% Last Modified by GUIDE v2.5 26-Jun-2018 16:26:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WaveMaker_OpeningFcn, ...
    'gui_OutputFcn',  @WaveMaker_OutputFcn, ...
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


% --- Executes just before WaveMaker is made visible.
function WaveMaker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WaveMaker (see VARARGIN)

% Choose default command line output for WaveMaker
handles.output = hObject;

handles.stSine=0.2;
handles.spSine=0.2;
handles.stCosine=0;
handles.spCosine=0;
handles.stNoise=0;
handles.spNoise=0;
handles.stNwaves=1;
handles.spNwaves=1;

saveops.fig2=0;
saveops.fig1=0;
saveops.mat=0;
saveops.vtk=1;

setappdata(0,'saveoptions',saveops);

try
    d = gpuDevice;
    set(handles.GPU_Acceleration,'Enable','on')
catch
    set(handles.GPU_Acceleration,'Enable','off')
end

set(handles.togglebuttonContinueSimulation,'Enable','off')

handles.LoadInitialConditionFileName='';
handles.LoadInitialConditionPathName='';
handles.UseExistingInitialCondition=0;
handles.InitialConditionCheck='New IC';

handles.New_LstBxValue = 1;

handles.ContinueSimulation = 0;


%Set all parameters in Dimensional, Nusselt and Skadov scaling
Reynolds = 40.8;
Kapitza = 3923;
Ct = 0;
k0x= 0.0729;
k0z= 0.0729;

set(handles.Reynolds,'String',sprintf('%.3g',Reynolds));
set(handles.Kapitza,'String',sprintf('%.3g',Kapitza));
set(handles.Ct,'String',sprintf('%.3g',Ct));

set(handles.k0x,'String',sprintf('%.3g',k0x));
set(handles.k0z,'String',sprintf('%.3g',k0z));



%Set dimensional values in Reynolds scaling
Shkadov_eta=(3*Reynolds)^(4/9)/Kapitza^(2/3);
Shkadov_delta=(3*Reynolds)^(11/9)/Kapitza^(1/3);
Shkadov_zeta=Ct*(3*Reynolds)^(2/9)/Kapitza^(1/3);
kappa=3*Reynolds/Shkadov_delta;
Shkadov_k0x=kappa*k0x;
Shkadov_k0z=kappa*k0z;

set(handles.Value_Shkadov_delta,'String',sprintf('%.3g',Shkadov_delta));
set(handles.Value_Shkadov_eta,'String',sprintf('%.3g',Shkadov_eta));
set(handles.Value_Shkadov_zeta,'String',sprintf('%.3g',Shkadov_zeta));
set(handles.Value_Shkadov_k0x,'String',sprintf('%.3g',Shkadov_k0x));
set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',Shkadov_k0z));
 %set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',Shkadov_k0z));
%Set dimensional Values

%Read Fixed Dimensional values
Dim_nu=1e-6;
Dim_gAbs=9.81;
Dim_rho=1000;

set(handles.Nu,'String',sprintf('%.3g',Dim_nu));
set(handles.gAbs,'String',sprintf('%.3g',Dim_gAbs));
set(handles.FluidProperties_rho,'String',sprintf('%.3g',Dim_rho));

if Ct>=0
    Dim_InclinationAngle=acot(Ct)*360/2/pi;
else
    Dim_InclinationAngle=180+acot(Ct)*360/2/pi;
end

Dim_gx=Dim_gAbs*sin(Dim_InclinationAngle/360*2*pi());

FluidProp_sigma=Kapitza*(Dim_rho*(Dim_gx)^(1/3)*Dim_nu^(4/3));
Dim_Filmthickness=(3*Reynolds*Dim_nu^2/Dim_gx)^(1/3);

FluidProp_Nu=str2num(get(handles.Nu,'String'));

lv=(FluidProp_Nu^2/Dim_gx)^(1/3);
hnDach=(3*Reynolds)^(1/3)*lv;

Dim_Lx=hnDach*2*pi/k0x;
Dim_Lz=hnDach*2*pi/k0z;

set(handles.FluidProperties_Sigma,'String',sprintf('%.3g',FluidProp_sigma));
set(handles.Dimensional_Filmthickness,'String',sprintf('%.3g',Dim_Filmthickness));
set(handles.Dimensional_InclinationAngle,'String',sprintf('%.3g',Dim_InclinationAngle));
set(handles.Lx,'String',sprintf('%.3g',Dim_Lx));
set(handles.Lz,'String',sprintf('%.3g',Dim_Lz));

%% determine most amplified wavelength in Nusselt scaling

Re=Reynolds;

Weber=Kapitza/((3*Re)^(2/3));

m=1;
for k=0.01:0.001:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersion relation of the Full second order S. 197/98 Falling Films (Kalliadasis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% omega^4 terms
A=-(12/715)*Re^3;

% omega^3 terms
B=-(90/143)*1i*Re^2+(108/5005)*Re^3*k-(2027/80080)*1i*Re^2*k^2;

% omega^2 terms
C=(54/13)*Re+(98/143)*1i*Re^2*k+(3231/3640)*Re*k^2+(27/5005)*Ct*Re^2*k^2-(612/65065)*Re^3*k^2+(3439/145600)*1i*Re^2*k^3+(27/5005)*Weber*Re^2*k^4;

% omega^1 terms
D= 3*1i-(522/143)*Re*k+(27/5)*1i*k^2+(12/65)*1i*Ct*Re*k^2-(26424/117117)*1i*Re^2*k^2-(2441/4004)*Re*k^3-(16/5005)*Ct*Re^2*k^3+(1104/715715)*Re^3*k^3-(4591/650650)*1i*Re^2*k^4+(12/65)*1i*Weber*Re*k^4-(16/5005)*Weber*Re^2*k^5;

% omega^0 terms
E=-3*1i*k+(498/715)*Re*k^2-Ct*k^2+(1368/65065)*1i*Re^2*k^3-(304/5005)*1i*Ct*Re*k^3-(12/5)*1i*k^3-Weber*k^4-(48/715715)*Re^3*k^4+(148/325325)*Ct*Re^2*k^4+(30993/320320)*Re*k^4-(304/5005)*1i*Weber*Re*k^5+(1773/2602600)*1i*Re^2*k^5+(148/325325)*Weber*Re^2*k^6;

p=[A B C D E];
    
x_full=[];
x_full=roots(p);
full_real(:,m)=real(x_full);
x_full=imag(x_full);
x_full_max(m)=max(x_full);
m=m+1;
end

pos_of_max_growth=find(x_full_max==max(x_full_max));
wavenumber_of_max_growth=pos_of_max_growth*0.001+0.01;
handles.most_amplified=wavenumber_of_max_growth;
guidata(hObject, handles);
set(handles.edit_most_amplified,'string',num2str(handles.most_amplified));



% Find all static text UICONTROLS whose 'Tag' starts with latex_

lbls = findobj(hObject,'-regexp','tag','latex_*');

for i=1:length(lbls)

      l = lbls(i);

      % Get current text, position and tag

      set(l,'units','normalized');

      s = strcat('\textit{$',get(l,'string'),'$}');

      
      p = get(l,'position');

      t = get(l,'tag');

      localparent=get(l,'Parent');
      
      % Remove the UICONTROL

      delete(l);

      % Replace it with a TEXT object 
      
      handles.laxis = axes('parent',localparent,'units','normalized','position',p,'visible','off');
      handles.(t) = text(0.5,0.5,s,'interpreter','latex','units','normalized','fontsize',12,'HorizontalAlignment','center');

end

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = WaveMaker_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SaveParameter.
function SaveParameter_Callback(hObject, eventdata, handles)

try
    %GET ParameterFileName
    a=get(handles.ParameterFileName,'String');
    parameter.Name=a;
    
    parameter.Model=get(handles.popupmenu4,'Value')
    parameter.Scaling=get(handles.Scaling,'Value')
    handles.New_LstBxValue=handles.Scaling;
    
    parameter.Reynolds=str2num(get(handles.Reynolds,'String'));
    parameter.Kapitza=str2num(get(handles.Kapitza,'String'));
    parameter.Ct=str2num(get(handles.Ct,'String'));
    parameter.h0=str2num(get(handles.h0,'String'));
    parameter.Nusselt_k0x=str2num(get(handles.k0x,'String'));
    parameter.Nusselt_k0z=str2num(get(handles.k0z,'String'));
    
    parameter.Shkadov_delta=str2num(get(handles.Value_Shkadov_delta,'String'));
    parameter.Shkadov_eta=str2num(get(handles.Value_Shkadov_eta,'String'));
    parameter.Shkadov_zeta=str2num(get(handles.Value_Shkadov_zeta,'String'));
    parameter.Shkadov_kappa=3*parameter.Reynolds/parameter.Shkadov_delta;
    parameter.Shkadov_k0x=str2num(get(handles.Value_Shkadov_k0x,'String'));
    parameter.Shkadov_k0z=str2num(get(handles.Value_Shkadov_k0z,'String'));   
    
    parameter.DimVar_nu=str2num(get(handles.Nu,'String'));
    parameter.DimVar_rho=str2num(get(handles.FluidProperties_rho,'String'));
    parameter.DimVar_sigma=str2num(get(handles.FluidProperties_Sigma,'String'));
    parameter.DimVar_gAbs=str2num(get(handles.gAbs,'String'));
    parameter.DimVar_hN=str2num(get(handles.Dimensional_Filmthickness,'String'));
    parameter.DimVar_Theta=str2num(get(handles.Dimensional_InclinationAngle,'String'));
    parameter.DimVar_Lx=str2num(get(handles.Lx,'String'));
    parameter.DimVar_Lz=str2num(get(handles.Lz,'String'));
     
    parameter.NNst=(get(handles.NNstPop,'Value'));
    parameter.NNsp=(get(handles.NNspPop,'Value'))-1;
    
    parameter.TMax=str2num(get(handles.TMax,'String'));
    parameter.PrtStep=str2num(get(handles.PrtStep,'String'));
    parameter.AlaisingRadius=str2num(get(handles.AlaisingRadius,'String'));
    parameter.stSine=handles.stSine;
    parameter.spSine=handles.spSine;
    parameter.stCosine=handles.stCosine;
    parameter.spCosine=handles.spCosine;
    parameter.stNoise=handles.stNoise;
    parameter.spNoise=handles.spNoise;
    parameter.stNwaves=handles.stNwaves;
    parameter.spNwaves=handles.spNwaves;
    
    parameter.GPU_Acceleration=get(handles.GPU_Acceleration,'Value');
    
    parameter.LoadInitialConditionFileName=handles.LoadInitialConditionFileName;
    parameter.LoadInitialConditionPathName=handles.LoadInitialConditionPathName;
    parameter.UseExistingInitialCondition=handles.UseExistingInitialCondition;
    
    handles.saveoptions=getappdata(0,'saveoptions');
    
    parameter.saveoptions.fig2=handles.saveoptions.fig2;
    parameter.saveoptions.fig1=handles.saveoptions.fig1;
    parameter.saveoptions.mat=handles.saveoptions.mat;
    parameter.saveoptions.vtk=handles.saveoptions.vtk;
    
    if ~exist('ParameterFiles', 'dir')
        mkdir('ParameterFiles');
    end
    
    save(['ParameterFiles/' parameter.Name], 'parameter')
    
    msgbox('Parameter Saving Completed','Parameter Loading')
catch
    msgbox('Parameter Saving Failed','Parameter Loading')
end

% --- Executes on button press in loadParameter.
function loadParameter_Callback(hObject, eventdata, handles)
try
    [FileName,PathName] = uigetfile('*.mat','Select the MATLAB code file','./ParameterFiles/');
    
    LoadedParameter=load([PathName FileName]);
    
    parameter=LoadedParameter.parameter;
    
    handles.New_LstBxValue=parameter.Scaling;
    set(handles.Scaling,'Value',handles.New_LstBxValue);
    
    set(handles.popupmenu4,'Value',parameter.Model);
    
    set(handles.Reynolds,'String',sprintf('%.3g',parameter.Reynolds));
    set(handles.Kapitza,'String',sprintf('%.3g',parameter.Kapitza));
    set(handles.Ct,'String',sprintf('%.3g',parameter.Ct));
    set(handles.h0,'String',sprintf('%.3g',parameter.h0));
    set(handles.k0x,'String',sprintf('%.3g',parameter.Nusselt_k0x));
    set(handles.k0z,'String',sprintf('%.3g',parameter.Nusselt_k0z));
    
    [~,name,~] = fileparts(FileName);
    set(handles.ParameterFileName,'String',name);
    
    set(handles.Value_Shkadov_delta,'String',sprintf('%.3g',parameter.Shkadov_delta));
    set(handles.Value_Shkadov_eta,'String',sprintf('%.3g',parameter.Shkadov_eta));
    set(handles.Value_Shkadov_zeta,'String',sprintf('%.3g',parameter.Shkadov_zeta));
    
    set(handles.GPU_Acceleration,'Value',parameter.GPU_Acceleration);
    
    set(handles.Value_Shkadov_k0x,'String',sprintf('%.3g',parameter.Shkadov_k0x));
    set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',parameter.Shkadov_k0z));
    
    
    set(handles.Nu,'String',sprintf('%.3g',parameter.DimVar_nu));
    set(handles.FluidProperties_rho,'String',sprintf('%.3g',parameter.DimVar_rho));
    set(handles.FluidProperties_Sigma,'String',sprintf('%.3g',parameter.DimVar_sigma));
    set(handles.gAbs,'String',sprintf('%.3g',parameter.DimVar_gAbs));
    set(handles.Dimensional_Filmthickness,'String',sprintf('%.3g',parameter.DimVar_hN));
    set(handles.Dimensional_InclinationAngle,'String',sprintf('%.3g',parameter.DimVar_Theta));
    set(handles.Lx,'String',sprintf('%.3g',parameter.DimVar_Lx));
    set(handles.Lz,'String',sprintf('%.3g',parameter.DimVar_Lz));
    
    set(handles.NNstPop,'Value',(parameter.NNst));
    set(handles.NNspPop,'Value',(parameter.NNsp+1));
    
    set(handles.TMax,'String',sprintf('%.3g',parameter.TMax));
    set(handles.PrtStep,'String',sprintf('%.3g',parameter.PrtStep));
    set(handles.AlaisingRadius,'String',sprintf('%.3g',parameter.AlaisingRadius));
    handles.stSine=parameter.stSine;
    handles.spSine=parameter.spSine;
    handles.stCosine=parameter.stCosine;
    handles.spCosine=parameter.spCosine;
    handles.stNoise=parameter.stNoise;
    handles.spNoise=parameter.spNoise;
    handles.stNwaves=parameter.stNwaves;
    handles.spNwaves=parameter.spNwaves;
    
    saveoptions.fig2=parameter.saveoptions.fig2;
    saveoptions.fig1=parameter.saveoptions.fig1;
    saveoptions.mat=parameter.saveoptions.mat;
    saveoptions.vtk=parameter.saveoptions.vtk;
    
    setappdata(0,'saveoptions',saveoptions);
    
    
    if parameter.UseExistingInitialCondition==true
        handles.InitialConditionCheck='Load IC';
    else
        handles.InitialConditionCheck='New IC';
    end
    
    %% determine most amplified wavelength in Nusselt scaling

Re=parameter.Reynolds;
Kapitza=parameter.Kapitza;
Ct=parameter.Ct;
Weber=Kapitza/((3*Re)^(2/3));

m=1;
for k=0.01:0.001:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersion relation of the Full second order S. 197/98 Falling Films (Kalliadasis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% omega^4 terms
A=-(12/715)*Re^3;

% omega^3 terms
B=-(90/143)*1i*Re^2+(108/5005)*Re^3*k-(2027/80080)*1i*Re^2*k^2;

% omega^2 terms
C=(54/13)*Re+(98/143)*1i*Re^2*k+(3231/3640)*Re*k^2+(27/5005)*Ct*Re^2*k^2-(612/65065)*Re^3*k^2+(3439/145600)*1i*Re^2*k^3+(27/5005)*Weber*Re^2*k^4;

% omega^1 terms
D= 3*1i-(522/143)*Re*k+(27/5)*1i*k^2+(12/65)*1i*Ct*Re*k^2-(26424/117117)*1i*Re^2*k^2-(2441/4004)*Re*k^3-(16/5005)*Ct*Re^2*k^3+(1104/715715)*Re^3*k^3-(4591/650650)*1i*Re^2*k^4+(12/65)*1i*Weber*Re*k^4-(16/5005)*Weber*Re^2*k^5;

% omega^0 terms
E=-3*1i*k+(498/715)*Re*k^2-Ct*k^2+(1368/65065)*1i*Re^2*k^3-(304/5005)*1i*Ct*Re*k^3-(12/5)*1i*k^3-Weber*k^4-(48/715715)*Re^3*k^4+(148/325325)*Ct*Re^2*k^4+(30993/320320)*Re*k^4-(304/5005)*1i*Weber*Re*k^5+(1773/2602600)*1i*Re^2*k^5+(148/325325)*Weber*Re^2*k^6;

p=[A B C D E];
    
x_full=[];
x_full=roots(p);
full_real(:,m)=real(x_full);
x_full=imag(x_full);
x_full_max(m)=max(x_full);
m=m+1;
end

pos_of_max_growth=find(x_full_max==max(x_full_max));
wavenumber_of_max_growth=pos_of_max_growth*0.001+0.01;
handles.most_amplified=wavenumber_of_max_growth;
guidata(hObject, handles);
set(handles.edit_most_amplified,'string',num2str(handles.most_amplified));

    clear a
    msgbox('Parameter Loading Completed','Parameter Loading')
catch
    msgbox('Parameter Loading Failed','Parameter Loading')
end

guidata(hObject, handles); % updates handles structure

function GPU_Acceleration_Callback(hObject, eventdata, handles)

function ParameterFileName_Callback(hObject, eventdata, handles)
% hObject    handle to ParameterFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ParameterFileName as text
%        str2double(get(hObject,'String')) returns contents of ParameterFileName as a double


% --- Executes during object creation, after setting all properties.
function ParameterFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParameterFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Reynolds_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Reynolds_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Kapitza_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function Kapitza_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ct_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Ct_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h0_Callback(hObject, eventdata, handles)
function h0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NNst_CreateFcn(hObject, eventdata, handles)
function NNst_Callback(hObject, eventdata, handles)

function NNsp_CreateFcn(hObject, eventdata, handles)
function NNsp_Callback(hObject, eventdata, handles)


function k0x_Callback(hObject, eventdata, handles)
% hObject    handle to k0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k0x as text
%        str2double(get(hObject,'String')) returns contents of k0x as a double




% --- Executes during object creation, after setting all properties.
function k0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k0z_Callback(hObject, eventdata, handles)
% hObject    handle to k0z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k0z as text
%        str2double(get(hObject,'String')) returns contents of k0z as a double


% --- Executes during object creation, after setting all properties.
function k0z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k0z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TMax_Callback(hObject, eventdata, handles)
% hObject    handle to TMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TMax as text
%        str2double(get(hObject,'String')) returns contents of TMax as a double


% --- Executes during object creation, after setting all properties.
function TMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PrtStep_Callback(hObject, eventdata, handles)
% hObject    handle to PrtStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PrtStep as text
%        str2double(get(hObject,'String')) returns contents of PrtStep as a double


% --- Executes during object creation, after setting all properties.
function PrtStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PrtStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CaseName_Callback(hObject, eventdata, handles)
% hObject    handle to CaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CaseName as text
%        str2double(get(hObject,'String')) returns contents of CaseName as a double


% --- Executes during object creation, after setting all properties.
function CaseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AlaisingRadius_Callback(hObject, eventdata, handles)
% hObject    handle to AlaisingRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AlaisingRadius as text
%        str2double(get(hObject,'String')) returns contents of AlaisingRadius as a double


% --- Executes during object creation, after setting all properties.
function AlaisingRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlaisingRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in PlotInitialC.
function PlotInitialC_Callback(hObject, eventdata, handles)

addpath(genpath('./scr/'));
if handles.UseExistingInitialCondition==true
    fileName=[handles.LoadInitialConditionPathName '/' handles.LoadInitialConditionFileName];
    data=load(fileName);
    U=data.Usave;

    N1=data.N1;
    N2=data.N2;
    M1=data.M1;
    M2=data.M2;
    
[datah,~,~]   = reconstruct_U(U,N1,N2,M1,M2);

    hmax=real(max(max(datah)));
    figure
    if(M2>2)
        surf([real(datah)],'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        view([-35 55]);
        axis([0 M2 0 M1 0 ceil(hmax)]);
    end
    
    if(M2==2)
        plot(0:1/(M1-1):1,flipud((real(datah(:,1)))),'k','LineWidth',2);
        view([90 0]);
        axis([0 1 0 hmax*1.1]);
        hold on
        h0=str2num(get(handles.h0,'String'));
        
        parameter.Shkadov_k0x=str2num(get(handles.Value_Shkadov_k0x,'String'));
        parameter.Shkadov_k0z=str2num(get(handles.Value_Shkadov_k0z,'String'));   
        
        k0x = parameter.Shkadov_k0x;
        k0z = parameter.Shkadov_k0z;
        
        parameter.AlaisingRadius=str2num(get(handles.AlaisingRadius,'String'));
        
        [kx,kz,~]=CPU_AliasingFilter(parameter.AlaisingRadius,N1,N2,M1,M2,k0x,k0z);
        
        
        [velocity,~,~,~,~,~,y,~,~,~]= reconstruct_velocity_fields(U(1,:),M1,M2,N1,N2,h0,kx,kz,1,1);
        
        psi=cumsum(squeeze(velocity(:,1,:))'-max(max(max(velocity))));
        max(max(max(velocity)));
        contour(0:1/(M1-1):1,y,-(fliplr(psi)),20);
    end
    
end

if  handles.UseExistingInitialCondition==false
    
    k0x=str2num(get(handles.k0x,'String'));
    k0z=str2num(get(handles.k0z,'String'));
    
    NN1=(get(handles.NNstPop,'Value'));
    NN1P=NN1+1;
    NN2=(get(handles.NNspPop,'Value'))-1;
    NN2P=NN2+1;
    
    stSine=handles.stSine;
    spSine=handles.spSine;
    stCosine=handles.stCosine;
    spCosine=handles.spCosine;
    stNoise=handles.stNoise
    spNoise=handles.spNoise
    stNwaves=handles.stNwaves;
    spNwaves=handles.spNwaves;
    
    h0=str2num(get(handles.h0,'String'));
    
    if(k0x~=0) Lx=2.*pi/k0x; end
    if(k0z~=0) Lz=2.*pi/k0z; end
    if(NN1==0) k0x=0;Lx=0;stSine=0;stCosine=0;NN1P=1; end
    if(NN2==0) k0z=0;Lz=0;spSine=0;spCosine=0;NN2P=1; end
    if(NN1==0||NN2==0) dim=2; else dim=3; end
    
    [datah,~,~,~] = Initial_condition(stNoise, Lz, Lx, h0, stSine, k0x, stCosine, stNwaves, k0z, spSine, spNwaves, spCosine, spNoise, NN1, NN2);
    
    N1=2^NN1; M1=2*N1; N2=2^NN2; M2=2*N2; N=4*N2*(N1+1);
    datah   = ifft2(datah*(M1*M2));
    figure
    surf([real(datah)],'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    view([-35 55]);
    axis([0 M2 0 M1 0 ceil(hmax)]);
end
%set(gcf, 'Toolbar', 'figure')

% Update handles structure
guidata(hObject, handles);

function Nu_Callback(hObject, eventdata, handles)
% hObject    handle to Nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nu as text
%        str2double(get(hObject,'String')) returns contents of Nu as a double


% --- Executes during object creation, after setting all properties.
function Nu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gAbs_Callback(hObject, eventdata, handles)
% hObject    handle to gAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gAbs as text
%        str2double(get(hObject,'String')) returns contents of gAbs as a double


% --- Executes during object creation, after setting all properties.
function gAbs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lx_Callback(hObject, eventdata, handles)
% hObject    handle to Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lx as text
%        str2double(get(hObject,'String')) returns contents of Lx as a double


% --- Executes during object creation, after setting all properties.
function Lx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lz_Callback(hObject, eventdata, handles)
% hObject    handle to Lz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lz as text
%        str2double(get(hObject,'String')) returns contents of Lz as a double


% --- Executes during object creation, after setting all properties.
function Lz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Model.
function Model_Callback(hObject, eventdata, handles)
% hObject    handle to Model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Model


% --- Executes during object creation, after setting all properties.
function Model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in RunSimplifiedModel.
% function RunSimplifiedModel_Callback(hObject, eventdata, handles)
% 
% runSimulationSimplifiedModel(handles)
% % hObject    handle to RunSimplifiedModel (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 

% --- Executes on button press in InputShkadov.
function InputShkadov_Callback(hObject, eventdata, handles)

defaultanswer={get(handles.Value_Shkadov_delta,'String'),get(handles.Value_Shkadov_eta,'String'),get(handles.Value_Shkadov_zeta,'String'),...
    get(handles.Value_Shkadov_k0x,'String'),get(handles.Value_Shkadov_k0z,'String')};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
name='Input for Peaks function';

[x]= inputdlg({'Reduced Reynolds number \delta','Viscous Dispersion number \eta','Reduced inclination number \zeta','k_{x}','k_{z}'},...
    'Input values in Shkadov scaling', 1,defaultanswer,options);

if exist('x')
Shkadov_delta = str2num(x{1});
Shkadov_eta = str2num(x{2});
Shkadov_zeta = str2num(x{3});
Shkadov_k0x=str2num(x{4});
Shkadov_k0z=str2num(x{5});

set(handles.Value_Shkadov_delta,'String',sprintf('%.3g',Shkadov_delta));
set(handles.Value_Shkadov_eta,'String',sprintf('%.3g',Shkadov_eta));
set(handles.Value_Shkadov_zeta,'String',sprintf('%.3g',Shkadov_zeta));
set(handles.Value_Shkadov_k0x,'String',sprintf('%.3g',Shkadov_k0x));
set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',Shkadov_k0z));

Reynolds=Shkadov_delta/(3* Shkadov_eta^(1/2));
Ct=Shkadov_zeta/(Shkadov_eta^(1/2));
Kapitza=Shkadov_delta^(2/3)/(Shkadov_eta^(11/6));


set(handles.Reynolds,'String',sprintf('%.3g',Reynolds));
set(handles.Kapitza,'String',sprintf('%.3g',Kapitza));
set(handles.Ct,'String',sprintf('%.3g',Ct));

kappa=3*Reynolds/Shkadov_delta;
set(handles.k0x,'String',sprintf('%.3g',Shkadov_k0x/kappa));
set(handles.k0z,'String',sprintf('%.3g',Shkadov_k0z/kappa));

%Set dimensional Values

%Read Fixed Dimensional values
Dim_nu=str2num(get(handles.Nu,'String'));
Dim_gAbs=str2num(get(handles.gAbs,'String'));
Dim_rho=str2num(get(handles.FluidProperties_rho,'String'));

if Ct>=0
    Dim_InclinationAngle=acot(Ct)*360/2/pi;
else
    Dim_InclinationAngle=180+acot(Ct)*360/2/pi;
end
Dim_gx=Dim_gAbs*sin(Dim_InclinationAngle/360*2*pi());

FluidProp_sigma=Kapitza*(Dim_rho*(Dim_gx)^(1/3)*Dim_nu^(4/3));
Dim_Filmthickness=(3*Reynolds*Dim_nu^2/Dim_gx)^(1/3);


FluidProp_Nu=str2num(get(handles.Nu,'String'));
Dim_gx=str2num(get(handles.gAbs,'String'));

lv=(FluidProp_Nu^2/Dim_gx)^(1/3);
hnDach=(3*Reynolds)^(1/3)*lv;

Dim_Lx=hnDach*2*pi/(Shkadov_k0x/kappa);
Dim_Lz=hnDach*2*pi/(Shkadov_k0z/kappa);


set(handles.FluidProperties_Sigma,'String',sprintf('%.3g',FluidProp_sigma));
set(handles.gAbs,'String',sprintf('%.3g',Dim_gAbs));
set(handles.Dimensional_Filmthickness,'String',sprintf('%.3g',Dim_Filmthickness));
set(handles.Dimensional_InclinationAngle,'String',sprintf('%.3g',Dim_InclinationAngle));
set(handles.Lx,'String',sprintf('%.3g',Dim_Lx));
set(handles.Lz,'String',sprintf('%.3g',Dim_Lz));

%% determine most amplified wavelength in Nusselt scaling

Re=Reynolds;

Weber=Kapitza/((3*Re)^(2/3));

m=1;
for k=0.01:0.001:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersion relation of the Full second order S. 197/98 Falling Films (Kalliadasis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% omega^4 terms
A=-(12/715)*Re^3;

% omega^3 terms
B=-(90/143)*1i*Re^2+(108/5005)*Re^3*k-(2027/80080)*1i*Re^2*k^2;

% omega^2 terms
C=(54/13)*Re+(98/143)*1i*Re^2*k+(3231/3640)*Re*k^2+(27/5005)*Ct*Re^2*k^2-(612/65065)*Re^3*k^2+(3439/145600)*1i*Re^2*k^3+(27/5005)*Weber*Re^2*k^4;

% omega^1 terms
D= 3*1i-(522/143)*Re*k+(27/5)*1i*k^2+(12/65)*1i*Ct*Re*k^2-(26424/117117)*1i*Re^2*k^2-(2441/4004)*Re*k^3-(16/5005)*Ct*Re^2*k^3+(1104/715715)*Re^3*k^3-(4591/650650)*1i*Re^2*k^4+(12/65)*1i*Weber*Re*k^4-(16/5005)*Weber*Re^2*k^5;

% omega^0 terms
E=-3*1i*k+(498/715)*Re*k^2-Ct*k^2+(1368/65065)*1i*Re^2*k^3-(304/5005)*1i*Ct*Re*k^3-(12/5)*1i*k^3-Weber*k^4-(48/715715)*Re^3*k^4+(148/325325)*Ct*Re^2*k^4+(30993/320320)*Re*k^4-(304/5005)*1i*Weber*Re*k^5+(1773/2602600)*1i*Re^2*k^5+(148/325325)*Weber*Re^2*k^6;

p=[A B C D E];
    
x_full=[];
x_full=roots(p);
full_real(:,m)=real(x_full);
x_full=imag(x_full);
x_full_max(m)=max(x_full);
m=m+1;
end

pos_of_max_growth=find(x_full_max==max(x_full_max));
wavenumber_of_max_growth=pos_of_max_growth*0.001+0.01;
handles.most_amplified=wavenumber_of_max_growth;
guidata(hObject, handles);
set(handles.edit_most_amplified,'string',num2str(handles.most_amplified));

end

% --- Executes on button press on input Nusselt scaling.
function InputReynolds_Callback(hObject, eventdata, handles)

defaultanswer={get(handles.Reynolds,'String'),get(handles.Kapitza,'String'),get(handles.Ct,'String'),...
    get(handles.k0x,'String'),get(handles.k0z,'String')};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

[x]= inputdlg({'Reynolds number Re','Kapitza number Ka','Inclination number Ct',...
    'k_{x}','k_{z}'},...
    'Input values in Nusselt scaling', 1,defaultanswer,options);

if exist('x')
Reynolds = str2num(x{1});
Kapitza = str2num(x{2});
Ct = str2num(x{3});
k0x=str2num(x{4});
k0z=str2num(x{5});

set(handles.Reynolds,'String',sprintf('%.3g',Reynolds));
set(handles.Kapitza,'String',sprintf('%.3g',Kapitza));
set(handles.Ct,'String',sprintf('%.3g',Ct));

set(handles.k0x,'String',sprintf('%.3g',k0x));
set(handles.k0z,'String',sprintf('%.3g',k0z));



%Set dimensional values in Reynolds scaling
Shkadov_eta=(3*Reynolds)^(4/9)/Kapitza^(2/3);
Shkadov_delta=(3*Reynolds)^(11/9)/Kapitza^(1/3);
Shkadov_zeta=Ct*(3*Reynolds)^(2/9)/Kapitza^(1/3);
kappa=3*Reynolds/Shkadov_delta;
Shkadov_k0x=kappa*k0x;
Shkadov_k0z=kappa*k0z;

set(handles.Value_Shkadov_delta,'String',sprintf('%.3g',Shkadov_delta));
set(handles.Value_Shkadov_eta,'String',sprintf('%.3g',Shkadov_eta));
set(handles.Value_Shkadov_zeta,'String',sprintf('%.3g',Shkadov_zeta));
set(handles.Value_Shkadov_k0x,'String',sprintf('%.3g',Shkadov_k0x));
set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',Shkadov_k0z));
 %set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',Shkadov_k0z));
%Set dimensional Values

%Read Fixed Dimensional values
Dim_nu=str2num(get(handles.Nu,'String'));
Dim_gAbs=str2num(get(handles.gAbs,'String'));
Dim_rho=str2num(get(handles.FluidProperties_rho,'String'));

if Ct>=0
    Dim_InclinationAngle=acot(Ct)*360/2/pi;
else
    Dim_InclinationAngle=180+acot(Ct)*360/2/pi;
end

Dim_gx=Dim_gAbs*sin(Dim_InclinationAngle/360*2*pi());

FluidProp_sigma=Kapitza*(Dim_rho*(Dim_gx)^(1/3)*Dim_nu^(4/3));
Dim_Filmthickness=(3*Reynolds*Dim_nu^2/Dim_gx)^(1/3);

FluidProp_Nu=str2num(get(handles.Nu,'String'));

lv=(FluidProp_Nu^2/Dim_gx)^(1/3);
hnDach=(3*Reynolds)^(1/3)*lv;

Dim_Lx=hnDach*2*pi/k0x;
Dim_Lz=hnDach*2*pi/k0z;

set(handles.FluidProperties_Sigma,'String',sprintf('%.3g',FluidProp_sigma));
set(handles.Dimensional_Filmthickness,'String',sprintf('%.3g',Dim_Filmthickness));
set(handles.Dimensional_InclinationAngle,'String',sprintf('%.3g',Dim_InclinationAngle));
set(handles.Lx,'String',sprintf('%.3g',Dim_Lx));
set(handles.Lz,'String',sprintf('%.3g',Dim_Lz));

%% determine most amplified wavelength in Nusselt scaling

Re=Reynolds;

Weber=Kapitza/((3*Re)^(2/3));

m=1;
for k=0.01:0.001:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersion relation of the Full second order S. 197/98 Falling Films (Kalliadasis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% omega^4 terms
A=-(12/715)*Re^3;

% omega^3 terms
B=-(90/143)*1i*Re^2+(108/5005)*Re^3*k-(2027/80080)*1i*Re^2*k^2;

% omega^2 terms
C=(54/13)*Re+(98/143)*1i*Re^2*k+(3231/3640)*Re*k^2+(27/5005)*Ct*Re^2*k^2-(612/65065)*Re^3*k^2+(3439/145600)*1i*Re^2*k^3+(27/5005)*Weber*Re^2*k^4;

% omega^1 terms
D= 3*1i-(522/143)*Re*k+(27/5)*1i*k^2+(12/65)*1i*Ct*Re*k^2-(26424/117117)*1i*Re^2*k^2-(2441/4004)*Re*k^3-(16/5005)*Ct*Re^2*k^3+(1104/715715)*Re^3*k^3-(4591/650650)*1i*Re^2*k^4+(12/65)*1i*Weber*Re*k^4-(16/5005)*Weber*Re^2*k^5;

% omega^0 terms
E=-3*1i*k+(498/715)*Re*k^2-Ct*k^2+(1368/65065)*1i*Re^2*k^3-(304/5005)*1i*Ct*Re*k^3-(12/5)*1i*k^3-Weber*k^4-(48/715715)*Re^3*k^4+(148/325325)*Ct*Re^2*k^4+(30993/320320)*Re*k^4-(304/5005)*1i*Weber*Re*k^5+(1773/2602600)*1i*Re^2*k^5+(148/325325)*Weber*Re^2*k^6;

p=[A B C D E];
    
x_full=[];
x_full=roots(p);
full_real(:,m)=real(x_full);
x_full=imag(x_full);
x_full_max(m)=max(x_full);
m=m+1;
end

pos_of_max_growth=find(x_full_max==max(x_full_max));
wavenumber_of_max_growth=pos_of_max_growth*0.001+0.01;
handles.most_amplified=wavenumber_of_max_growth;
guidata(hObject, handles);
set(handles.edit_most_amplified,'string',num2str(handles.most_amplified));

end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)


defaultanswer={get(handles.Nu,'String'),get(handles.FluidProperties_rho,'String'),get(handles.FluidProperties_Sigma,'String')};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

[x]= inputdlg({'Kinematic viscosity \nu','Density \rho','Surface tension \sigma'},...
    'Insert fluid properties', 1,defaultanswer,options);

if exist('x')
FluidProp_Nu = str2num(x{1});
FluidProp_rho = str2num(x{2});
FluidProp_sigma = str2num(x{3});

set(handles.Nu,'String',sprintf('%.3g',FluidProp_Nu));
set(handles.FluidProperties_rho,'String',sprintf('%.3g',FluidProp_rho));
set(handles.FluidProperties_Sigma,'String',sprintf('%.3g',FluidProp_sigma));
end


% --- Executes on button press in SetDimensionalProperties.
function SetDimensionalProperties_Callback(hObject, eventdata, handles)


defaultanswer={get(handles.gAbs,'String'),get(handles.Dimensional_Filmthickness,'String'),get(handles.Dimensional_InclinationAngle,'String'),...
    get(handles.Lx,'String'),get(handles.Lz,'String')};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

[x]= inputdlg({'Gravitational acceleration | g |','Nusselt flat film thickness h_N','Inclination angle \theta','Streamwise wavelength L_x','Spanwise wavelength L_z'},...
    'Insert dimensional parameters', 1,defaultanswer,options);

if exist('x')
Dim_gAbs = str2num(x{1});
Dim_Filmthickness = str2num(x{2});
Dim_InclinationAngle = str2num(x{3});
Dim_Lx=str2num(x{4});
Dim_Lz=str2num(x{5});

set(handles.gAbs,'String',sprintf('%.3g',Dim_gAbs));
set(handles.Dimensional_Filmthickness,'String',sprintf('%.3g',Dim_Filmthickness));
set(handles.Dimensional_InclinationAngle,'String',sprintf('%.3g',Dim_InclinationAngle));
set(handles.Lx,'String',sprintf('%.3g',Dim_Lx));
set(handles.Lz,'String',sprintf('%.3g',Dim_Lz));

%Change Nusselt Scaling Values
FluidProp_Nu=str2num(get(handles.Nu,'String'));
FluidProp_rho=str2num(get(handles.FluidProperties_rho,'String'));
FluidProp_sigma=str2num(get(handles.FluidProperties_Sigma,'String'));
sin(Dim_InclinationAngle*2*pi/360);

Dim_gx=Dim_gAbs*sin(Dim_InclinationAngle/360*2*pi());

Reynolds=Dim_gx*Dim_Filmthickness^3/(3*FluidProp_Nu^2);
Kapitza=FluidProp_sigma/(FluidProp_rho*(Dim_gx)^(1/3)*FluidProp_Nu^(4/3));

Ct = cot(Dim_InclinationAngle/360*2*pi);

lv=(FluidProp_Nu^2/Dim_gx)^(1/3);
hnDach=(3*Reynolds)^(1/3)*lv;

k0x=2*pi/(Dim_Lx/hnDach); %In Nusselt scaling
k0z=2*pi/(Dim_Lz/hnDach); %In Nusselt scaling

set(handles.Reynolds,'String',sprintf('%.3g',Reynolds));
set(handles.Kapitza,'String',sprintf('%.3g',Kapitza));
set(handles.Ct,'String',sprintf('%.3g',Ct));
set(handles.k0x,'String',sprintf('%.3g',k0x));
set(handles.k0z,'String',sprintf('%.3g',k0z));

%Change Shkadov Scaling Values
Shkadov_eta=(3*Reynolds)^(4/9)/Kapitza^(2/3);
Shkadov_delta=(3*Reynolds)^(11/9)/Kapitza^(1/3);
Shkadov_zeta=Ct*(3*Reynolds)^(2/9)/Kapitza^(1/3);

set(handles.Value_Shkadov_delta,'String',sprintf('%.3g',Shkadov_delta));
set(handles.Value_Shkadov_eta,'String',sprintf('%.3g',Shkadov_eta));
set(handles.Value_Shkadov_zeta,'String',sprintf('%.3g',Shkadov_zeta));
kappa=3*Reynolds/Shkadov_delta;
set(handles.Value_Shkadov_k0x,'String',sprintf('%.3g',k0x*kappa));
set(handles.Value_Shkadov_k0z,'String',sprintf('%.3g',k0z*kappa));

parameter.Reynolds=str2num(get(handles.Reynolds,'String'));
parameter.Kapitza=str2num(get(handles.Kapitza,'String'));
parameter.Ct=str2num(get(handles.Ct,'String'));

%% determine most amplified wavelength in Nusselt scaling

Re=parameter.Reynolds;
Ct=parameter.Ct;
Kapitza=parameter.Kapitza;

Weber=parameter.Kapitza/((3*Re)^(2/3));

m=1;
for k=0.01:0.001:4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersion relation of the Full second order S. 197/98 Falling Films (Kalliadasis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% omega^4 terms
A=-(12/715)*Re^3;

% omega^3 terms
B=-(90/143)*1i*Re^2+(108/5005)*Re^3*k-(2027/80080)*1i*Re^2*k^2;

% omega^2 terms
C=(54/13)*Re+(98/143)*1i*Re^2*k+(3231/3640)*Re*k^2+(27/5005)*Ct*Re^2*k^2-(612/65065)*Re^3*k^2+(3439/145600)*1i*Re^2*k^3+(27/5005)*Weber*Re^2*k^4;

% omega^1 terms
D= 3*1i-(522/143)*Re*k+(27/5)*1i*k^2+(12/65)*1i*Ct*Re*k^2-(26424/117117)*1i*Re^2*k^2-(2441/4004)*Re*k^3-(16/5005)*Ct*Re^2*k^3+(1104/715715)*Re^3*k^3-(4591/650650)*1i*Re^2*k^4+(12/65)*1i*Weber*Re*k^4-(16/5005)*Weber*Re^2*k^5;

% omega^0 terms
E=-3*1i*k+(498/715)*Re*k^2-Ct*k^2+(1368/65065)*1i*Re^2*k^3-(304/5005)*1i*Ct*Re*k^3-(12/5)*1i*k^3-Weber*k^4-(48/715715)*Re^3*k^4+(148/325325)*Ct*Re^2*k^4+(30993/320320)*Re*k^4-(304/5005)*1i*Weber*Re*k^5+(1773/2602600)*1i*Re^2*k^5+(148/325325)*Weber*Re^2*k^6;

p=[A B C D E];
    
x_full=[];
x_full=roots(p);
full_real(:,m)=real(x_full);
x_full=imag(x_full);
x_full_max(m)=max(x_full);
m=m+1;
end

pos_of_max_growth=find(x_full_max==max(x_full_max));
wavenumber_of_max_growth=pos_of_max_growth*0.001+0.01;
handles.most_amplified=wavenumber_of_max_growth;
guidata(hObject, handles);
set(handles.edit_most_amplified,'string',num2str(handles.most_amplified));

end



% --- Executes during object creation, after setting all properties.
function Value_Shkadov_delta_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Value_Shkadov_eta_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Value_Shkadov_eta_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Value_Shkadov_zeta_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Value_Shkadov_zeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Shkadov_zeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Value_Shkadov_k0x_Callback(hObject, eventdata, handles)
% hObject    handle to Value_Shkadov_k0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_Shkadov_k0x as text
%        str2double(get(hObject,'String')) returns contents of Value_Shkadov_k0x as a double


% --- Executes during object creation, after setting all properties.
function Value_Shkadov_k0x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Shkadov_k0x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Value_Shkadov_k0z_Callback(hObject, eventdata, handles)
% hObject    handle to Value_Shkadov_k0z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_Shkadov_k0z as text
%        str2double(get(hObject,'String')) returns contents of Value_Shkadov_k0z as a double


% --- Executes during object creation, after setting all properties.
function Value_Shkadov_k0z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Shkadov_k0z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FluidProperties_rho_Callback(hObject, eventdata, handles)
% hObject    handle to FluidProperties_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FluidProperties_rho as text
%        str2double(get(hObject,'String')) returns contents of FluidProperties_rho as a double


% --- Executes during object creation, after setting all properties.
function FluidProperties_rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FluidProperties_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FluidProperties_Sigma_Callback(hObject, eventdata, handles)
% hObject    handle to FluidProperties_Sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FluidProperties_Sigma as text
%        str2double(get(hObject,'String')) returns contents of FluidProperties_Sigma as a double


% --- Executes during object creation, after setting all properties.
function FluidProperties_Sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FluidProperties_Sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dimensional_Filmthickness_Callback(hObject, eventdata, handles)
% hObject    handle to Dimensional_Filmthickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dimensional_Filmthickness as text
%        str2double(get(hObject,'String')) returns contents of Dimensional_Filmthickness as a double


% --- Executes during object creation, after setting all properties.
function Dimensional_Filmthickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dimensional_Filmthickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dimensional_InclinationAngle_Callback(hObject, eventdata, handles)
% hObject    handle to Dimensional_InclinationAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dimensional_InclinationAngle as text
%        str2double(get(hObject,'String')) returns contents of Dimensional_InclinationAngle as a double


% --- Executes during object creation, after setting all properties.
function Dimensional_InclinationAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dimensional_InclinationAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetInitialConditions.
function SetInitialConditions_Callback(hObject, eventdata, handles)
% hObject    handle to SetInitialConditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%

defaultanswer={sprintf('%.3g',handles.stSine), sprintf('%.3g',handles.stCosine), sprintf('%.3g',handles.spSine),...
    sprintf('%.3g',handles.spCosine), sprintf('%.3g',handles.stNoise), sprintf('%.3g',handles.spNoise), sprintf('%.3g',handles.stNwaves), sprintf('%.3g',handles.spNwaves)};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

[x]= inputdlg({'Streamwise amplitude sine','Streamwise amplitude cosine',...
    'Spanwise amplitude sine','Spanwise amplitude cosine',...
    '2-D noise','Spanwise noise',...
    'Number of streamwise waves','Number of spanwise waves'},...
    'Input values in dimensional numbers', 1,defaultanswer,options);

if exist('x')
handles.stSine=str2num(x{1});
handles.stCosine=str2num(x{2});
handles.spSine=str2num(x{3});
handles.spCosine=str2num(x{4});
handles.stNoise=str2num(x{5});
handles.spNoise=str2num(x{6});
handles.stNwaves=str2num(x{7});
handles.spNwaves=str2num(x{8});

handles.UseExistingInitialCondition=false;
handles.InitialConditionCheck='New IC';

guidata(hObject, handles); % updates handles structure
end



% --- Executes on button press in LoadInitialCondition.
function LoadInitialCondition_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*_U_Full2ndOrder.mat','Select the MATLAB code file','./');

try
handles.LoadInitialConditionFileName=FileName;
handles.LoadInitialConditionPathName=PathName;
handles.UseExistingInitialCondition=true;
handles.InitialConditionCheck='Load IC';

    fileName=[handles.LoadInitialConditionPathName '/' handles.LoadInitialConditionFileName];
    data=load(fileName);

    N1=data.N1;
    N2=data.N2;
    
    set(handles.NNstPop,'Value',log2(N1));
    set(handles.NNspPop,'Value',log2(N2)+1);

guidata(hObject, handles); % updates handles structure
catch
    msgbox('Loading initial condition failed','Loading initial condition')
end
% hObject    handle to LoadInitialCondition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function NNstPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NNstPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function NNspPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NNspPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Scaling.
function Scaling_Callback(hObject, eventdata, handles)
% hObject    handle to Scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Scaling contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Scaling

% get parameter data
parameter.Reynolds=str2num(get(handles.Reynolds,'String'));
parameter.Kapitza=str2num(get(handles.Kapitza,'String'));
parameter.Ct=str2num(get(handles.Ct,'String'));
parameter.h0=str2num(get(handles.h0,'String'));
parameter.Nusselt_k0x=str2num(get(handles.k0x,'String'));
parameter.Nusselt_k0z=str2num(get(handles.k0z,'String'));

parameter.Shkadov_delta=str2num(get(handles.Value_Shkadov_delta,'String'));
parameter.Shkadov_eta=str2num(get(handles.Value_Shkadov_eta,'String'));
parameter.Shkadov_zeta=str2num(get(handles.Value_Shkadov_zeta,'String'));
parameter.Shkadov_kappa=3*parameter.Reynolds/parameter.Shkadov_delta;
parameter.Shkadov_k0x=str2num(get(handles.Value_Shkadov_k0x,'String'));
parameter.Shkadov_k0z=str2num(get(handles.Value_Shkadov_k0z,'String'));

parameter.DimVar_nu=str2num(get(handles.Nu,'String'));
parameter.DimVar_rho=str2num(get(handles.FluidProperties_rho,'String'));
parameter.DimVar_sigma=str2num(get(handles.FluidProperties_Sigma,'String'));
parameter.DimVar_gAbs=str2num(get(handles.gAbs,'String'));
parameter.DimVar_hN=str2num(get(handles.Dimensional_Filmthickness,'String'));
parameter.DimVar_Theta=str2num(get(handles.Dimensional_InclinationAngle,'String'));
parameter.DimVar_Lx=str2num(get(handles.Lx,'String'));
parameter.DimVar_Lz=str2num(get(handles.Lz,'String'));

% determine dependend variables

DimVar_gx=parameter.DimVar_gAbs*sin(parameter.DimVar_Theta/360*2*pi());
lv=(parameter.DimVar_nu^2/DimVar_gx)^(1/3);
tv=(parameter.DimVar_nu/(DimVar_gx^2))^(1/3);
hnDach=(3*parameter.Reynolds)^(1/3)*lv;

% get old and new Scaling

handles.Old_LstBxValue = handles.New_LstBxValue;
handles.New_LstBxValue=get(handles.Scaling,'Value');

% change display of Tend and Tprint according to old and new scaling

Tend_old=str2num(get(handles.TMax,'String'));
PrtStep_old=str2num(get(handles.PrtStep,'String'));

if handles.Old_LstBxValue == 1 && handles.New_LstBxValue == 2
    set(handles.TMax,'String',sprintf('%.3g',Tend_old*hnDach/(tv*lv)));
    set(handles.PrtStep,'String',sprintf('%.3g',PrtStep_old*hnDach/(tv*lv)));
end

if handles.Old_LstBxValue == 1 && handles.New_LstBxValue == 3
    set(handles.TMax,'String',sprintf('%.3g',Tend_old*hnDach/(tv*lv*parameter.Shkadov_kappa)));
    set(handles.PrtStep,'String',sprintf('%.3g',PrtStep_old*hnDach/(tv*lv*parameter.Shkadov_kappa)));
end

if handles.Old_LstBxValue == 2 && handles.New_LstBxValue == 1
    set(handles.TMax,'String',sprintf('%.3g',Tend_old*(tv*lv)/hnDach));
    set(handles.PrtStep,'String',sprintf('%.3g',PrtStep_old*(tv*lv)/hnDach));
end

if handles.Old_LstBxValue == 2 && handles.New_LstBxValue == 3
    set(handles.TMax,'String',sprintf('%.3g',Tend_old/parameter.Shkadov_kappa));
    set(handles.PrtStep,'String',sprintf('%.3g',PrtStep_old/parameter.Shkadov_kappa));
end

if handles.Old_LstBxValue == 3 && handles.New_LstBxValue == 1
    set(handles.TMax,'String',sprintf('%.3g',Tend_old*tv*lv/hnDach*parameter.Shkadov_kappa));
    set(handles.PrtStep,'String',sprintf('%.3g',PrtStep_old*tv*lv/hnDach*parameter.Shkadov_kappa));
end

if handles.Old_LstBxValue == 3 && handles.New_LstBxValue == 2
    set(handles.TMax,'String',sprintf('%.3g',Tend_old*parameter.Shkadov_kappa));
    set(handles.PrtStep,'String',sprintf('%.3g',PrtStep_old*parameter.Shkadov_kappa));
end

guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function Scaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

% if(get(handles.togglebuttonContinueSimulation,'BackgroundColor'))~= 'red'

set(handles.togglebutton1,'string','Stop','BackgroundColor','red');
set(handles.togglebuttonContinueSimulation,'Enable','off')
set(handles.togglebuttonContinueSimulation,'BackgroundColor',[0.94 0.94 0.94]);
set(handles.togglebuttonContinueSimulation,'ForegroundColor','black');

pause(.001);

a=get(handles.ParameterFileName,'String');
parameter.Name=a;

parameter.Model=get(handles.popupmenu4,'Value');
parameter.Scaling=get(handles.Scaling,'Value');

parameter.Reynolds=str2num(get(handles.Reynolds,'String'));
parameter.Kapitza=str2num(get(handles.Kapitza,'String'));
parameter.Ct=str2num(get(handles.Ct,'String'));
parameter.h0=str2num(get(handles.h0,'String'));
parameter.Nusselt_k0x=str2num(get(handles.k0x,'String'));
parameter.Nusselt_k0z=str2num(get(handles.k0z,'String'));

parameter.Shkadov_delta=str2num(get(handles.Value_Shkadov_delta,'String'));
parameter.Shkadov_eta=str2num(get(handles.Value_Shkadov_eta,'String'));
parameter.Shkadov_zeta=str2num(get(handles.Value_Shkadov_zeta,'String'));
parameter.Shkadov_kappa=3*parameter.Reynolds/parameter.Shkadov_delta;
parameter.Shkadov_k0x=str2num(get(handles.Value_Shkadov_k0x,'String'));
parameter.Shkadov_k0z=str2num(get(handles.Value_Shkadov_k0z,'String'));

parameter.DimVar_nu=str2num(get(handles.Nu,'String'));
parameter.DimVar_rho=str2num(get(handles.FluidProperties_rho,'String'));
parameter.DimVar_sigma=str2num(get(handles.FluidProperties_Sigma,'String'));
parameter.DimVar_gAbs=str2num(get(handles.gAbs,'String'));
parameter.DimVar_hN=str2num(get(handles.Dimensional_Filmthickness,'String'));
parameter.DimVar_Theta=str2num(get(handles.Dimensional_InclinationAngle,'String'));
parameter.DimVar_Lx=str2num(get(handles.Lx,'String'));
parameter.DimVar_Lz=str2num(get(handles.Lz,'String'));

parameter.NNst=(get(handles.NNstPop,'Value'));
parameter.NNsp=(get(handles.NNspPop,'Value'))-1;

parameter.TMax=str2num(get(handles.TMax,'String'));
parameter.PrtStep=str2num(get(handles.PrtStep,'String'));
parameter.AlaisingRadius=str2num(get(handles.AlaisingRadius,'String'));
parameter.stSine=handles.stSine;
parameter.spSine=handles.spSine;
parameter.stCosine=handles.stCosine;
parameter.spCosine=handles.spCosine;
parameter.stNoise=handles.stNoise;
parameter.spNoise=handles.spNoise;
parameter.stNwaves=handles.stNwaves;
parameter.spNwaves=handles.spNwaves;

parameter.GPU_Acceleration=get(handles.GPU_Acceleration,'Value');

parameter.LoadInitialConditionFileName=handles.LoadInitialConditionFileName;
parameter.LoadInitialConditionPathName=handles.LoadInitialConditionPathName;
parameter.UseExistingInitialCondition=handles.UseExistingInitialCondition;

addpath(genpath('./scr/'));
parameter.continue=0;

parameter.azFig1 = -35;
parameter.elFig1 = 55;
[parameter.UEnd, parameter.TEnd]=runSimulation(handles,hObject,parameter);
save(['./' parameter.Name '/ParameterFile'],'parameter');
set(handles.togglebuttonContinueSimulation,'Enable','on')
set(handles.togglebutton1,'string','Start','BackgroundColor','blue');
% end
% --- Executes on button press in SaveOptions.
function SaveOptions_Callback(hObject, eventdata, handles)
% hObject    handle to SaveOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create figure

f=figure('Name','Save Options','NumberTitle','off','units','pixels','position',[200,200,200,100],...
    'toolbar','none','menu','none');
movegui(f,'center');
% Create save options checkboxes
c1=uicontrol('style','checkbox','units','pixels',...
    'position',[30,40,60,15],'string','fig1.fig','tag','fig1');
c2=uicontrol('style','checkbox','units','pixels',...
    'position',[120,40,60,15],'string','fig2.mat','tag','fig2');
c3=uicontrol('style','checkbox','units','pixels',...
    'position',[30,70,60,15],'string','fig1.mat','tag','mat');
c4=uicontrol('style','checkbox','units','pixels',...
    'position',[120,70,60,15],'string','.vtk','tag','vtk');
% Create OK pushbutton
uicontrol('style','pushbutton','units','pixels',...
    'position',[65,5,70,20],'string','OK',...
    'callback',@yourCallback);

handles.saveoptions=getappdata(0,'saveoptions');

set(c1,'Value',handles.saveoptions.fig1)    ;
set(c2,'Value',handles.saveoptions.fig2)    ;
set(c3,'Value',handles.saveoptions.mat)    ;
set(c4,'Value',handles.saveoptions.vtk)    ;


function yourCallback (hObject, eventdata, handles)
localVars=struct();
localVars.vtk=hObject.Parent.Children(2).Value;
localVars.mat=hObject.Parent.Children(3).Value;
localVars.fig2=hObject.Parent.Children(4).Value;
localVars.fig1=hObject.Parent.Children(5).Value;

setappdata(0,'saveoptions',localVars);

close


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CreateStruct.Interpreter = 'latex';
CreateStruct.WindowStyle = 'modal';
CreateStruct.FontSize = 15;
%h=msgbox('Z = X^2 + Y^2','Value',CreateStruct);

h=msgbox({'\underline{Nusselt scaling:}'...
    ,'','$Re = \frac{g \sin \theta h_N^3}{3 \nu^2}$','',...
    '$Ka = \frac{\sigma}{\rho (g \sin \theta)^{1/3} \nu^{4/3}}$','',...
    '$Ct = \cot(\theta)$','',...
    '$k_x = \frac{2 \pi}{L_x/h_N}$','',...
    '$k_z = \frac{2 \pi}{L_z/h_N}$',''},'Nusselt',CreateStruct);
pos=get(h,'Position');
set(h, 'position',[pos(1),pos(2),150,230]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 13 ); %makes text bigger

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Scaling.
function Scaling_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Enable', 'Off');
handles.Old_LstBxValue=get(handles.Scaling,'Value');
guidata(hObject, handles);


% --- Executes on button press in SetDimensionalProperties.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to SetDimensionalProperties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CreateStruct.Interpreter = 'latex';
CreateStruct.WindowStyle = 'modal';
CreateStruct.FontSize = 15;
%h=msgbox('Z = X^2 + Y^2','Value',CreateStruct);

h=msgbox({'\underline{Shkadov scaling:}'...
    ,'','$\delta = \frac{(3 Re)^{11/9}}{Ka^{1/3}}$','',...
    '$\eta = \frac{(3 Re)^{4/9}}{Ka^{2/3}}$','',...
    '$\zeta = \frac{Ct(3 Re)^{2/9}}{Ka^{1/3}}$',''...
    '$k_x =  \frac{k_{x,Nu}}{\sqrt{\eta}}$',''...
    '$k_z =  \frac{k_{z,Nu}}{\sqrt{\eta}}$',''...
    },'Shkadov',CreateStruct);
pos=get(h,'Position');
set(h, 'position',[pos(1),pos(2),150,230]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 13 ); %makes text bigger


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CreateStruct.Interpreter = 'latex';
CreateStruct.WindowStyle = 'modal';
CreateStruct.FontSize = 13;
%h=msgbox('Z = X^2 + Y^2','Value',CreateStruct);

h=msgbox({'\underline{Save options}'...
    ,'','\textbf{fig1.mat}: save film height data in chosen scaling','',...
    '\textbf{.vtk}: save velocity and alpha field compatible to paraview','',...
    '\textbf{fig1.fig}: save figure 1 as matlab figure file','',...
    '\textbf{fig2.mat}: save data shown in figure 2 as table',''},'Save options',CreateStruct);
pos=get(h,'Position');
set(h, 'position',[pos(1),pos(2),320,180]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 11 ); %makes text bigger



% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CreateStruct.Interpreter = 'latex';
CreateStruct.WindowStyle = 'modal';
CreateStruct.FontSize = 11;
%h=msgbox('Z = X^2 + Y^2','Value',CreateStruct);

h=msgbox({'$h_0$ relates the Reynolds number based on the Nusselt film thickness'...
'(periodic domain) to the Reynolds number based on the actual flow rate in an open domain.','',...
'A value around 0.9 gives a good relation between the two Reynolds numbers in many cases.',''},'Parameter definition',CreateStruct);
pos=get(h,'Position');
set(h, 'position',[pos(1),pos(2),350,150]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 10 ); %makes text bigger





% --- Executes on button press in pbExit.
function pbExit_Callback(hObject, eventdata, handles)
% hObject    handle to pbExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)


% --- Executes on selection change in NNspPop.
function NNspPop_Callback(hObject, eventdata, handles)
% hObject    handle to NNspPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NNspPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NNspPop


% --- Executes on selection change in NNstPop.
function NNstPop_Callback(hObject, eventdata, handles)
% hObject    handle to NNstPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NNstPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NNstPop


% --- Executes on button press in togglebuttonContinueSimulation.
function togglebuttonContinueSimulation_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonContinueSimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonContinueSimulation
set(handles.togglebuttonContinueSimulation,'string','Stop','BackgroundColor','red');
set(handles.togglebuttonContinueSimulation,'string','Stop','ForegroundColor','white');
set(handles.togglebutton1,'Enable','off');
a=get(handles.ParameterFileName,'String');
parameter.Name=a;
load(['./' parameter.Name '/ParameterFile'])
   [parameter.azFig1,parameter.elFig1] = view(handles.axes2);
    
parameter.continue = 1;
[parameter.UEnd, parameter.TEnd]=runSimulation(handles,hObject,parameter);
save(['./' parameter.Name '/ParameterFile'],'parameter');
set(handles.togglebuttonContinueSimulation,'string','Continue','BackgroundColor','blue');
set(handles.togglebutton1,'Enable','on');

% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit_most_amplified_Callback(hObject, eventdata, handles)
% hObject    handle to edit_most_amplified (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_most_amplified as text
%        str2double(get(hObject,'String')) returns contents of edit_most_amplified as a double


% --- Executes during object creation, after setting all properties.
function edit_most_amplified_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_most_amplified (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CreateStruct.Interpreter = 'latex';
CreateStruct.WindowStyle = 'modal';
CreateStruct.FontSize = 11;
%h=msgbox('Z = X^2 + Y^2','Value',CreateStruct);

h=msgbox({'$\mathrm{k_{x,max,Nu}}$ is the wavenumber in Nusselt scaling, which '...
'corresponds to the most amplified streamwise wavelength.','',...
'$\mathrm{k_{x,max,Nu}}$ is derived from the linear dispersion relation of the full Second-Order Model.',''},'Parameter definition',CreateStruct);
pos=get(h,'Position');
set(h, 'position',[pos(1),pos(2),280,120]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 10 ); %makes text bigger
