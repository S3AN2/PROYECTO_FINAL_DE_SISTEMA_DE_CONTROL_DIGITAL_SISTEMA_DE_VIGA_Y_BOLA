function varargout = controldigitalb(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @controldigitalb_OpeningFcn, ...
                   'gui_OutputFcn',  @controldigitalb_OutputFcn, ...
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


% --- Executes just before controldigitalb is made visible.
function controldigitalb_OpeningFcn(hObject, eventdata, handles, varargin)
axes(handles.axes4);
path = '../assets/FIEE.jpg';
img = imread(path);
imshow(img);
axis off;
axes(handles.axes3);
path = '../assets/UNMSM.jpg';
img = imread(path);
imshow(img);
axis off;
axes(handles.axes5);
path = '../assets/controlpid.png';
img = imread(path);
imshow(img);
axis off;


find_system('Name','../SIMULINK/PARTEB');
open_system('../SIMULINK/PARTEB');
handles.output = hObject;
guidata(hObject, handles);


% UIWAIT makes controldigitalb wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = controldigitalb_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global CC GG HH DD  
global ecua
global controlPID
global fundis
global ecuPID
global control
global controloca
global Gdiscreta
global GD
global CPID
global sys
stopt=get(handles.stopia,'string');
param_masa2=get(handles.masa,'string');
param_radio2=get(handles.radio,'string');
brazo2=get(handles.brazo,'string');
longitud2=get(handles.longitud,'string');
momento2=get(handles.momento,'string');
periodo2=get(handles.periodo,'string');
% tiempo2=get(handles.periodo,'string');
param_MP=get(handles.MP,'string');
param_TSS=get(handles.TSS,'string');
Mp1=eval(param_MP);
Mp=Mp1/100;
Tss=eval(param_TSS);

radio1=eval(param_radio2);
masa1=eval(param_masa2);
brazo1=eval(brazo2);
longitud1=eval(longitud2);
momento1=eval(momento2);
periodo1=eval(periodo2);
stopi=eval(stopt)
% tiempo1=eval(tiempo2);
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
Cf=C
sistema=ss(A,B,C,D)
[num,den]=ss2tf (A,B,C,D);
funcion=tf(num,den)
num3=num
den3=den
Hd=c2d(funcion,t);
funcion_de_transferencia_discreta=Hd;
GD=Hd;
num4=GD.num
num5=cell2mat(num4)

den4=GD.den
den5=cell2mat(den4)

[G,HD,Cd,Dd]=c2dm(A,B,C,D,t,'zoh');
G;
HD;
Cd;
Dd;
H=HD;
C=Cd;
D=Dd;
%PERIODO
T=0.01
%FACTOR DE AMORTIGUAMIENTO
ep = sqrt((log(Mp))^2/(pi^2+(log(Mp))^2)); %factor de amortiguamiento
%FRECUENCIA NATURAL
wn = 4/(ep*Tss);%frecuencia natural del sistema

%
%
%
%Ecuacion deseada
Gd = tf([wn^2],[1 2*wn*ep wn^2]);

%Calcular los polos del sistema deseado
Polos_deseados = pole(Gd) ;%Raices complejas 
Polo_3 = 10*real(Polos_deseados);
j=ep;
Wn=wn;
% %calculo de los parametros del PID
% k = num;
% Ps = den;
% %polinomio deseado
% p1 = conv([1 -real(Polos_deseados(1))],[1 -real(Polos_deseados(1))]) + [0 0 imag(Polos_deseados(1))^2];
% %ecuacion deseada
% Ecuacion_deseada= conv([1 -Polo_3(1)],p1);
% Ecuacion_deseada_del_PID=Ecuacion_deseada;
% %valores del PID
% Kc = (Ecuacion_deseada(3)-Ps(3))/k
% Ti = (k*Kc)/Ecuacion_deseada(4)
% Td = (Ecuacion_deseada(2)-Ps(2))/(k*Kc)
Ec1=[1 (2*j*Wn) (Wn^2)];
%Hallando las raíces
z=real(roots(Ec1));
%Matriz auxiliarpara hallar la parte real de los polos dominantes
a=[0 1];
w=10*(a*z);%Tercer polo es 5 veces 
Ec2=[1 -w];
%Ecuacion total
Ec_to=conv(Ec1,Ec2);
Ecuacion_deseada_del_PID=Ec_to;
KD=Ec_to(2)/0.21
KP=Ec_to(3)/0.21
KI=Ec_to(4)/0.21

N=char('KD : ','KP : ','KI : ');
% N=char('KC : ','Ti : ','Td : ');
% R=[Kc;Ti;Td];
R=[KD;KP;KI];
S=num2str(R);
PARAMETROS_DEL_CONTROLADOR_PID=[N,S];
    Ti=KP/KI;
    Td=KD/KP;
    Kpimc=KP;
    %
    s=tf('s');
    G1=funcion;
    sys1=parallel(Kpimc,Kpimc*1/(Ti*s));
    GPID=parallel(sys1,Kpimc*s*Td);
    G7=feedback(GPID*G1,1);
    CPID=c2d(G7,T,'zoh')
%
% poolos=[(Polos_deseados(1)) (Polos_deseados(2))];
%
[A,B,C,D]=tf2ss(num,den);
funcion_de_transferencia=tf(num,den);

%
%
%
% PM(Z)=1+p1Z^(-1)+p2Z^(-2)
p1=-2*exp(-ep*wn*T)*cos(wn*T*sqrt(1-ep^2))
p2=exp(-2*ep*wn*T)
%raices de los polos  PM(Z)=1+p1Z^(-1)+p2Z^(-2) 

ecdeseada=[1 p1 p2]
Ecuacion_deseada_del_control_por_localizacion=num2str(ecdeseada);
raices=roots(ecdeseada)

%POLOS A USAR 
J=[raices(1) raices(2)]

%MATRIZ DE GANANCIA METODO ACKERMAN
K=acker(G,H,J)

%ENTRADA
u=1
%SALIDA
Yss=C*inv(eye(2)-G+H*K)*H*K(1)*u
%METODO GRAFICO
GG=G-H*K;
HH=H*K(1)
%HH=[0;1300];
CC=[1 0];
DD=[0];
sys=ss(GG,HH,CC,DD,T)

%controlabilidad
Coo=ctrb(G,H);
dime=size(Coo);
Dimension=dime(1);
rangoC=rank(Coo);
Rango_de_la_controlabilidad=rangoC;

%
N1=char('Controlabilidad : ');
N3=char('Dimension: ');
N2=char('Rango: ');
S2=mat2str(Coo);
S4=mat2str(Dimension);
S3=mat2str(rangoC);
hhh=[N1,S2];
HHH=[N2,S3];
kkk=[N3,S4];
if Dimension==rangoC
   tex= ("El sistema es controlable")
else
    tex= ("No es controlable")
end
N4=tex;
Controlabilidad_del_sistema=strvcat(hhh,kkk,HHH,N4);

%SIMULINK
%
set_param('PARTEB','StopTime',num2str(stopi));
%PID
set_param('PARTEB/PID','P',num2str(KP));
set_param('PARTEB/PID','I',num2str(KI));
set_param('PARTEB/PID','D',num2str(KD));
%FUNCION DE TRANSFERENCIA
set_param('PARTEB/tf1','Numerator',mat2str(num5),'Denominator',mat2str(den5));
%
%CONTROLADOR POR LOCALIZACION DE POLOS
set_param('PARTEB/G4','Gain',mat2str(H));
set_param('PARTEB/G5','Gain',mat2str(G));
set_param('PARTEB/G6','Gain',mat2str(CC));
set_param('PARTEB/G','Gain',num2str(K(1)));
set_param('PARTEB/G1','Gain',num2str(K(2)));
%
%
%
%

%variables
 
 ecua=evalc('Ecuacion_deseada_del_control_por_localizacion');
 controlPID=evalc('PARAMETROS_DEL_CONTROLADOR_PID');

 fundis=evalc('funcion_de_transferencia_discreta');
 ecuPID=evalc('Ecuacion_deseada_del_PID');

 control=evalc('Controlabilidad_del_sistema');

 controloca=evalc('K');
 Gdiscreta=evalc('funcion_de_transferencia_discreta')







% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
close(controldigitalb)
controldigital


function MP_Callback(hObject, eventdata, handles)
% hObject    handle to MP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MP as text
%        str2double(get(hObject,'String')) returns contents of MP as a double


% --- Executes during object creation, after setting all properties.
function MP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TSS_Callback(hObject, eventdata, handles)
% hObject    handle to TSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TSS as text
%        str2double(get(hObject,'String')) returns contents of TSS as a double


% --- Executes during object creation, after setting all properties.
function TSS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
set_param(gcs,'SimulationCommand','Start');
open_system('PARTEB/Scope1');


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
set_param(gcs,'SimulationCommand','Start');
open_system('PARTEB/Scope2');


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
set_param(gcs,'SimulationCommand','Start');
open_system('PARTEB/Scope3');


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
% hObject    handle to periodo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of periodo as text
%        str2double(get(hObject,'String')) returns contents of periodo as a double


% --- Executes during object creation, after setting all properties.
function periodo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to periodo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)

global GD
global CPID
global sys

inf=get(hObject,'Value');
gos=get(hObject,'String');



switch inf
 case 2
set(handles.axes6); % Establece los ejes de graficación
axes(handles.axes6);
step(GD);grid on;
title('Función de transferencia discreta')
 case 3
set(handles.axes6); % Establece los ejes de graficación
axes(handles.axes6);
step(CPID);grid on;
title('CONTROLADOR PID')
  case 4
set(handles.axes6); % Establece los ejes de graficación
axes(handles.axes6);
step(sys);grid on;
title('CONTROLADOR POR LOCALIZACION DE POLOS')
  case 5
set(handles.axes6); % Establece los ejes de graficación
axes(handles.axes6);
step(CPID);
hold on
step(sys);grid on;
hold off
title('COMPARACIÓN AMBAS GRÁFICAS')


end
guidata(hObject,handles);

function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)


global fundis
global ecuPID
global ecua
global control
global controlPID
global controloca


inf=get(hObject,'Value');
gos=get(hObject,'String');
switch inf
 case 2
 set(handles.text122,'string',fundis)
 case 3
 set(handles.text122,'string',ecuPID)
  case 4
 set(handles.text122,'string',ecua)
  case 5
 set(handles.text122,'string',control)
  case 6
 set(handles.text122,'string',controlPID)
  case 7
 set(handles.text122,'string',controloca)
    
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stopia_Callback(hObject, eventdata, handles)
% hObject    handle to stopia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stopia as text
%        str2double(get(hObject,'String')) returns contents of stopia as a double


% --- Executes during object creation, after setting all properties.
function stopia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
