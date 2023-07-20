function varargout = controldigital(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @controldigital_OpeningFcn, ...
                   'gui_OutputFcn',  @controldigital_OutputFcn, ...
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


% --- Executes just before controldigital is made visible.
function controldigital_OpeningFcn(hObject, eventdata, handles, varargin)
axes(handles.axes1);
path = '../assets/FIEE.jpg';
img = imread(path);
imshow(img);
axis off;
axes(handles.axes2);
path = '../assets/UNMSM.jpg';
img = imread(path);
imshow(img);
axis off;
axes(handles.axes3);
path = '../assets/bola.jpeg';
img = imread(path);
imshow(img);
axis off;


find_system('Name','../SIMULINK/controol');
open_system('../SIMULINK/controol');

handles.output = hObject;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = controldigital_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

function slider1_Callback(hObject, eventdata, handles)




function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)




function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)




function listbox1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radiobutton2_Callback(~, eventdata, handles)

function botondigital_Callback(hObject, eventdata, handles)


% --- Executes on selection change in lista.
function lista_Callback(hObject, eventdata, handles)

global sistema2s
global funcions
global sistemas
global Hds
inf=get(hObject,'Value');
gos=get(hObject,'String');
switch inf
 case 1
 set(handles.resultado,'string',funcions)
 case 2
 set(handles.resultado,'string',sistemas)
  case 3
 set(handles.resultado,'string',Hds)
  case 4
 set(handles.resultado,"string",sistema2s)    
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function lista_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lista (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function masa_Callback(hObject, eventdata, handles)
% hObject    handle to masa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of masa as text
%        str2double(get(hObject,'String')) returns contents of masa as a double


% --- Executes during object creation, after setting all properties.
function masa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to masa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radio_Callback(hObject, eventdata, handles)
% hObject    handle to radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radio as text
%        str2double(get(hObject,'String')) returns contents of radio as a double


% --- Executes during object creation, after setting all properties.
function radio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in graficas.
function graficas_Callback(hObject, eventdata, handles)

global funcion
global sistema
global Hd

inf=get(hObject,'Value');



switch inf
 case 1
set(handles.axes9); % Establece los ejes de graficación
axes(handles.axes9);
step(funcion);grid on;
title('Función de transferencia')
 case 2
set(handles.axes9); % Establece los ejes de graficación
axes(handles.axes9);
step(sistema);grid on;
title('Espacio de Estados Continuo')
  case 3
set(handles.axes9); % Establece los ejes de graficación
axes(handles.axes9);
step(Hd,'--')
axis([0 50 0 300]);%limite de imagen
grid on;
title('Espacio de Estados Discreto')
  case 4
set(handles.axes9); % Establece los ejes de graficación
axes(handles.axes9);
step(Hd,'--',funcion,'g');
axis([0 50 0 300]);
grid on;
title('Espacio de Estados Discreto y Continuo')
  case 5
set(handles.axes9); % Establece los ejes de graficación
axes(handles.axes9);
%simulacion en el dominio del tiempo resuelta por espacio estado.
t2=0:1:50;
y=0.1047729*t2.^2;
plot(t2,y)
axis([0 50 0 300]);
grid on;
title('Dominio del tiempo resuelta por espacio estado')
  case 6
set(handles.axes9); % Establece los ejes de graficación
axes(handles.axes9);
%simulacion en el dominio del tiempo resuelta por espacio estado.
t2=0:1:50;
y=0.1047729*t2.^2;
hold on
plot(t2,y)

step(funcion)

hold off
grid on;
title('Dominio del tiempo y funcion de transferencia .')
 
end

function graficas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to graficas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function brazo_Callback(hObject, eventdata, handles)
% hObject    handle to brazo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of brazo as text
%        str2double(get(hObject,'String')) returns contents of brazo as a double


% --- Executes during object creation, after setting all properties.
function brazo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brazo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function longitud_Callback(hObject, eventdata, handles)
% hObject    handle to longitud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of longitud as text
%        str2double(get(hObject,'String')) returns contents of longitud as a double


% --- Executes during object creation, after setting all properties.
function longitud_CreateFcn(hObject, eventdata, handles)
% hObject    handle to longitud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function momento_Callback(hObject, eventdata, handles)
% hObject    handle to momento (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of momento as text
%        str2double(get(hObject,'String')) returns contents of momento as a double


% --- Executes during object creation, after setting all properties.
function momento_CreateFcn(hObject, eventdata, handles)
% hObject    handle to momento (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function periodo_Callback(hObject, eventdata, handles)
% hObject    handle to text12213 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text12213 as text
%        str2double(get(hObject,'String')) returns contents of text12213 as a double


% --- Executes during object creation, after setting all properties.
function periodo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text12213 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tiempo_Callback(hObject, eventdata, handles)
% hObject    handle to tiempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tiempo as text
%        str2double(get(hObject,'String')) returns contents of tiempo as a double


% --- Executes during object creation, after setting all properties.
function tiempo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tiempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ejecutar.
function ejecutar_Callback(hObject, eventdata, handles)
global sistema
global funcion
global sistema2
global Hd
global sistema2s
global funcions
global sistemas
global Hds
syms z
param_masa2=get(handles.masa,'string');
param_radio2=get(handles.radio,'string');
brazo2=get(handles.brazo,'string');
longitud2=get(handles.longitud,'string');
momento2=get(handles.momento,'string');
periodo2=get(handles.periodo,'string');
tiempo2=get(handles.tiempo,'string');

radio1=eval(param_radio2);
masa1=eval(param_masa2);
brazo1=eval(brazo2);
longitud1=eval(longitud2);
momento1=eval(momento2);
periodo1=eval(periodo2);
tiempo1=eval(tiempo2);


m=masa1;
g=-9.8;
d=brazo1;
L=longitud1;
J=momento1*10^(-6);
R=radio1;
t=periodo1;

A=[0 1;0 0];
B=[0;-(m*g*d)/(L*((J/R^2)+m))];
C=[1 0];
D=[0];
sistema=ss(A,B,C,D)
[num,den]=ss2tf (A,B,C,D);
funcion=tf(num,den)

Hd=c2d(funcion,t)
[G,HD,Cd,Dd]=c2dm(A,B,C,D,t,'zoh');
G;
HD;
Cd;
Dd;
sistema2=ss(G,HD,Cd,Dd);
sistema2s=evalc("sistema2");
funcions=evalc("funcion");
sistemas=evalc("sistema");
Hds=evalc("Hd");
set_param('controol/FT','Numerator',mat2str(num),'Denominator',mat2str(den));
set_param('controol/estadoC','A',mat2str(A),'B',mat2str(B),'C',mat2str(C),'D',mat2str(D));
set_param('controol/estadoD','A',mat2str(G),'B',mat2str(HD),'C',mat2str(Cd),'D',mat2str(Dd),"Sampletime",num2str(periodo1));
set_param('controol','StopTime',num2str(tiempo1));
set_param(gcs,'SimulationCommand','Start');
open_system('controol/Scope');
% --- Executes during object creation, after setting all properties.
function resultado_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in animacion.
function animacion_Callback(hObject, eventdata, handles)
ballbeam


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
close(controldigital)
controldigitalb
