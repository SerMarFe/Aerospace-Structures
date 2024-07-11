clear all; close all; clc;
%%
addpath("code")
addpath("structureSolver_adapted_from_codeP1\")
format long

%% 1) PREPROCESS
% 1.1 Input data (define your input parameters here)
input_data

data.ni = 3; % Degrees of freedom per node
data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom
data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);

% Material properties matrix
E1=150e9;
E2=78e9;
A1=pi*(1.25e-3)^2;
A2=pi*((7e-3)^2-(5e-3)^2);
density1=1500;
density2=2600;
Atela=3.75;
mTela = 1.55*Atela;
mPIC = 85;
Cl=2.6;
Cd=1.45;
g=9.81;
sigmaYbar=240;
sigmaYcable=150;

m = [% Each column corresponds to a material property (area, Young's modulus, etc.)
    E1    A1    0  density1
    E2    A2    0  density2
];

% Mass and Inertia
Mel = getMass(data,x,Tn,m,Tm);
Mtotal = sum(Mel) + mTela + mPIC;
mNodos = masaNodos(data,Mel,Tn);
mNodos_original = mNodos;
mNodos([3 4 5 6]) = mNodos([3 4 5 6]) + [mTela/3 mTela/3 mTela/6 mTela/6]';
mNodos([1 2]) = mNodos([1 2]) + mPIC/2;

posCM = (mNodos'*x)/Mtotal;
xCM = (x-posCM);

Ixx = mNodos'*(xCM(:,2).^2+xCM(:,3).^2);
Iyy = mNodos'*(xCM(:,1).^2+xCM(:,3).^2);
Izz = mNodos'*(xCM(:,1).^2+xCM(:,2).^2);
Ixy = mNodos'*(xCM(:,1).*xCM(:,2));
Ixz = mNodos'*(xCM(:,1).*xCM(:,3));
Iyz = mNodos'*(xCM(:,2).*xCM(:,3));

I = [Ixx -Ixy -Ixz
     -Ixy Iyy -Iyz
     -Ixz -Iyz Izz];

%% 1.2
% Tenemos un sistema en equilibrio excepto en z, lo cual implica suma de fuerzas (2 eq, en 'x','y') y de
% momentos (1 eq, y) igual a 0 (con la de momentos en 'x' y en 'z' hemos obtenido que las reacciones en 'x' y en 'z' son iguales).
% Incógnitas, 3 (Fpiloto en 1 y 2 son iguales para con la suma de momentos alrededor de z, x; y la velocidad Vx)
% En realidad Rx y Rz son funciones del tiempo ya que dependen de la
% velocidad en z que es función del tiempo. Para la resolución hemos
% sustituido Rx y Rz en función de vz para resolver la ecuación
% diferencial, hayando vz(t), a partir de la cual ya podemos conocer Rx(t)
% y Rz(t)

% Primero, encontrar Vx para saber las fuerzas de lift y drag

syms vz(t) Rx Rz

L = 1/2*1.225*Atela*Cl*vz^2;
D = 1/2*1.225*Atela*Cd*vz^2;

% Suma de fuerzas (x)
eq1=2*Rx-D==0;

% Suma de fuerzas (z)
eq2=2*Rz+L-Mtotal*g==Mtotal*diff(vz,t);

% Suma de momentos (y)
M12 = cross(x(1,:)-posCM,[2*Rx;0;2*Rz]);
% M12 = 2*Rz*abs(x(1,1)-posCM(1)) + 2*Rx*abs(x(1,3)-posCM(3));

M3 = cross(x(3,:)-posCM,[-D/3;0;L/3]);
% M3 = L/3*(x(3,1)-posCM(1)) + D/3*(x(3,3)-posCM(3));

M4 = cross(x(4,:)-posCM,[-D/3;0;L/3]);
% M4 = L/3*(x(4,1)-posCM(1)) + D/3*(x(4,3)-posCM(3));

M5 = cross(x(5,:)-posCM,[-D/6;0;L/6]);
% M5 = L/6*(x(5,1)-posCM(1)) + D/6*(x(5,3)-posCM(3));

M6 = cross(x(6,:)-posCM,[-D/6;0;L/6]);
% M6 = L/6*(x(6,1)-posCM(1)) + D/6*(x(6,3)-posCM(3));

eq3 = (M12 + M3 + M4 + M5 + M6)*[0;1;0] == 0;

%%

deq3_vz = subs(eq3,[Rx,Rz],[solve(eq1,Rx),solve(eq2,Rz)]);

vz_sol = dsolve(deq3_vz,vz(0)==0);
vz_func = matlabFunction(vz_sol);

Rx = subs(solve(eq1,Rx),vz,vz_sol);
Rx_func = matlabFunction(Rx);
Rz = subs(solve(eq2,Rz),[diff(vz,t),vz],[diff(vz_sol,t),vz_sol]);
Rz_func = matlabFunction(Rz);

L_func = matlabFunction(subs(L,vz,vz_sol));
D_func = matlabFunction(subs(D,vz,vz_sol));

%% 1.2 Gráficos
t = 0:0.01:3;
figure('Name','Reacciones del piloto');
plot(t,Rz_func(t),'-b', t,Rx_func(t),'m-');
title('Reacciones del piloto');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Fuerza de reacci\''on del piloto [N]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$R_z$','$R_x$','Interpreter','latex', 'Location','best'); %"Guany cr\'{i}tic "+sprintf("($K_u = %.4f$)",k_u);
grid minor, grid on;

%%
figure('Name','Velocidad de descenso');
plot(t,vz_func(t),'-c');
title('Velocidad de descenso');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Velocidad de descenso ($v_z$) [m/s]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$v_z$','Interpreter','latex', 'Location','best'); %"Guany cr\'{i}tic "+sprintf("($K_u = %.4f$)",k_u);
grid minor, grid on;


%% 1.3
% Ahora la fuerza de lift L se incrementa según una función del tiempo,
% solamente en el nodo 3, lo cual implica que el momento varía y la suma de
% fuerzas también. La reacción del piloto también se incrementa y se
% reparte en la misma forma entre los nodos 1 y 2. Por lo tanto, tenemos
% que la fuerza en x sigue teniendo la misma expresión, nada ha cambiado ni
% en Rx ni en el drag D.
% Puesto que la racha de viento se produce en t=3s la solución de la
% ecuación que obtendremos será aplicable solamente para t>=3s, para t<3
% usaremos la solución del apartado 1.2, donde no había la perturbación

syms t Rx Rz

vg = 2.4;
c1 = 4.512; c2 = 6.0996; c3 = 1.9218;

tg = 3; Delta_tg = 0.7; tr = tg + 0.5; Delta_tr = 2;
tf = tr+Delta_tr+2;

D = 1/2*1.225*Atela*Cd*vz^2;
L = 1/2*1.225*Atela*Cl*vz^2;
Lg = 1/2*1.225*Atela*Cl*vg^2;
Delta_Lg = Lg*sin(pi*(t-tg)/Delta_tg);

Delta_R = -Lg*(c3*(t-tr)^3 - c2*(t-tr)^2 + c1*(t-tr));

Delta_R_func = matlabFunction(Delta_R);
Delta_Lg_func = matlabFunction(Delta_Lg);

% Suma de fuerzas (x)
eq1=2*Rx-D==0;

% Suma de fuerzas (z)
eq2_1=Delta_Lg +2*Rz+L-Mtotal*g==Mtotal*diff(vz,t); % tg<=t<=tr
eq2_2=Delta_Lg+Delta_R +2*Rz+L-Mtotal*g==Mtotal*diff(vz,t); % tr<t<=tg+Delta_tg
eq2_3=Delta_R +2*Rz+L-Mtotal*g==Mtotal*diff(vz,t); % tg+Delta_tg<t<=tr+Delta_tr
eq2_4=2*Rz+L-Mtotal*g==Mtotal*diff(vz,t); % tr+Delta_tr<=t<tf

% Suma de momentos (y) (condición de equilibrio para vuelo estándard)
eq3 = (M12 + M3 + M4 + M5 + M6)*[0;1;0] == 0;


%% 1.3 Resolución numérica: Rx, Rz y vz
syms t v
dt = 0.01;
t_vec = tg:dt:tf;
Rx_vec = zeros(size(t_vec));
Rz_vec = zeros(size(t_vec));
vz_vec = zeros(size(t_vec));
vz_vec(1) = vz_func(tg);
Rx_vec(1) = Rx_func(tg);
Rz_vec(1) = Rz_func(tg);

for i=2:length(t_vec)
    if t_vec(i)<tr
        sols = solve(subs([eq1,eq2_1,eq3],[diff(vz,t),vz,t],[(v-vz_vec(i-1))/dt,v,t_vec(i)]));
    elseif t_vec(i)<tg+Delta_tg
        sols = solve(subs([eq1,eq2_2,eq3],[diff(vz,t),vz,t],[(v-vz_vec(i-1))/dt,v,t_vec(i)]));
    elseif t_vec(i)<tr+Delta_tr
        sols = solve(subs([eq1,eq2_3,eq3],[diff(vz,t),vz,t],[(v-vz_vec(i-1))/dt,v,t_vec(i)]));
    elseif t_vec(i)<tr+tf
        sols = solve(subs([eq1,eq2_4,eq3],[diff(vz,t),vz,t],[(v-vz_vec(i-1))/dt,v,t_vec(i)]));
    end

    v_ = sols.v;
    Rx_ = sols.Rx;
    Rz_ = sols.Rz;

    vz_vec(i) = double(v_(1));  % la solución buena (la segunda no tiene sentido);
    Rx_vec(i) = double(Rx_(1)); % la solución buena (la segunda no tiene sentido);
    Rz_vec(i) = double(Rz_(1)); % la solución buena (la segunda no tiene sentido);
end

%% Cálculo del lift
L_vec = double(subs(L,vz^2,vz_vec.^2));
D_vec = double(subs(D,vz^2,vz_vec.^2));


%% 1.4 Cálculo de la velocidad angular y aceleración angular
% Suma de momentos (y) para encontrar la aceleración angular
syms t omega_y(t) omega

M12 = cross(x(1)-posCM,[0;0;Delta_R])*[0;1;0];
M3 = cross(x(3)-posCM,[0;0;Delta_Lg])*[0;1;0];

eq4_1 = M3 == Iyy*diff(omega_y,t); % tg<=t<=tr
eq4_2 = M12 + M3 == Iyy*diff(omega_y,t); % tr<t<=tg+Delta_tg
eq4_3 = M12 == Iyy*diff(omega_y,t); % tg+Delta_tg<t<=tr+Delta_tr
eq4_4 = 0 == Iyy*diff(omega_y,t); % tr+Delta_tr<=t<tf

dt = 0.01;
t_vec = tg:dt:tf;
omega_y_vec = zeros(size(t_vec));
omega_y_vec(1) = 0;

for i=2:length(t_vec)
    if t_vec(i)<tr
        omega_y_vec(i) = solve(subs(eq4_1,[diff(omega_y,t),t],[(omega-omega_y_vec(i-1))/dt,t_vec(i)]));
    elseif t_vec(i)<tg+Delta_tg
        omega_y_vec(i) = solve(subs(eq4_2,[diff(omega_y,t),t],[(omega-omega_y_vec(i-1))/dt,t_vec(i)]));
    elseif t_vec(i)<tr+Delta_tr
        omega_y_vec(i) = solve(subs(eq4_3,[diff(omega_y,t),t],[(omega-omega_y_vec(i-1))/dt,t_vec(i)]));
    elseif t_vec(i)<tr+tf
        omega_y_vec(i) = solve(subs(eq4_4,[diff(omega_y,t),t],[(omega-omega_y_vec(i-1))/dt,t_vec(i)]));
    end
end


%% 1.4 Cálculo de la aceleración total de los nodos
t_vec1 = 0:dt:tg;
t_vec_full = [t_vec1 t_vec(2:end)];

vz_full = [vz_func(t_vec1) vz_vec(2:end)];
v_vec = [zeros(2,length(vz_full)); vz_full];
dvdt_vec = diff(v_vec,1,2)/dt;

omega_y_full = [zeros(1,length(t_vec1)) omega_y_vec(2:end)];
omega_vec = [zeros(1,length(omega_y_full)); omega_y_full; zeros(1,length(omega_y_full))];
domegadt_vec = diff(omega_vec,1,2)/dt;

rCM_vec = xCM; % Vector de posiciones de las masas de los nodos respecto el centro de masas (calculado anteriormente en el primer apartado);

a_vec = zeros(3,length(t_vec_full),data.nnod); % 3 filas (x,y,z), tantas columnas como tiempos guardados, tanta profundidad como número de nodos (la primera capa es la del primer nodo, etc)

for i=1:data.nnod
    for n=1:(length(t_vec_full)-1)
        a_vec(:,n,i) = dvdt_vec(:,n) + cross(domegadt_vec(:,n),rCM_vec(i,:))' + cross(omega_vec(:,n),cross(omega_vec(:,n),rCM_vec(i,:)))';
    end
end


%% 1.3 Gráficos: Lift
% Para el Lift
% tramo de 0 a tg: (no hay ráfaga)
t = 0:0.01:tg;
plot(t,L_func(t),'-g');
hold on;
% tramo de tg a tf
plot(t_vec,L_vec,'-g')

% Para el Drag
% tramo de 0 a tg: (no hay ráfaga)
t = 0:0.01:tg;
plot(t,D_func(t),'-r');
hold on;
% tramo de tg a tf
plot(t_vec,D_vec,'-r')

% Para la reacción del piloto
% Rz
% tramo de 0 a tr: (no hay reacción)
t = 0:0.01:tg;
plot(t,Rz_func(t),'-b');

% tramo de tg a tr: (no hay reacción, sí ráfaga)
t = tg:dt:tr;
plot(t,Rz_vec(t_vec<=tr),'-b');

% tramo de tr a tr+Delta_tr: (hay reacción)
t = tr:dt:tr+Delta_tr;
plot(t,Rz_vec(t_vec<=tr+Delta_tr & t_vec>=tr) + Delta_R_func(t),'-b');

% a partir de tr+Delta_tr: (no hay reacción)
t = tr+Delta_tr:0.01:tf;
plot(t,Rz_vec(t_vec>=tr+Delta_tr),'-b');

% Rx
% tramo de 0 a tg
t = 0:dt:tg;
plot(t,Rx_func(t),'-m')
% tramo de tg a tf
plot(t_vec,Rx_vec,'-m')
hold off;

title('Fuerzas en el ala delta');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Fuerzas [N]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$L$','','D','','$R_z$','','','','$R_x$','Interpreter','latex', 'Location','best'); %"Guany cr\'{i}tic "+sprintf("($K_u = %.4f$)",k_u);
grid minor, grid on;


%% 1.3 Gráficos: velocidad
% tramo de 0 a tg: (no hay ráfaga)
t = 0:0.01:tg;
figure('Name','Velocidad de descenso');
plot(t,vz_func(t),'-c');
hold on;

% tramo de tg a tr
t = tg:0.01:tf;
plot(t,vz_vec,'-c');
hold off;

title('Velocidad de descenso');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Velocidad de descenso ($v_z$) [m/s]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$v_z$','Interpreter','latex', 'Location','best');
grid minor, grid on;

%% 1.3 Gráficos: velocidad angular y aceleración angular
% tramo de 0 a tf:
figure('Name','Velocidad angular');
plot([0 t_vec],[0 omega_y_vec],'-g');

title('Velocidad angular en $y$','interpreter','latex');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Velocidad angular en $y$ ($\omega_y$) [rad/s]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$\omega_y$','Interpreter','latex', 'Location','best');
grid minor, grid on;
%%
% tramo de 0 a tf:
figure('Name','Aceleración angular');
plot(t_vec_full(1:end-1),domegadt_vec(2,:),'-b');
title('Aceleración angular en $y$','interpreter','latex');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Velocidad angular en $y$ ($\alpha_y$) [rad/s]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$\alpha_y$','Interpreter','latex', 'Location','best');
grid minor, grid on;

%% 1.4 Gráficos: aceleración total de cada nodo
figure('Name','Aceleración total de cada nodo');
plot(t_vec_full,reshape(a_vec(1,:,:),length(t_vec_full),data.nnod));legend('1','2','3','4','5','6')
title('Aceleraci\''on $x$ total de cada nodo','interpreter','latex');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Aceleraci\''on total en $x$ ($\omega_y$) [m/s$^2$]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

figure('Name','Aceleración total de cada nodo');
plot(t_vec_full,reshape(a_vec(2,:,:),length(t_vec_full),data.nnod));legend('1','2','3','4','5','6')
title('Aceleraci\''on $y$ total de cada nodo','interpreter','latex');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Aceleraci\''on total en $y$ ($\omega_y$) [m/s$^2$]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

figure('Name','Aceleración total de cada nodo');
plot(t_vec_full,reshape(a_vec(3,:,:),length(t_vec_full),data.nnod));legend('1','2','3','4','5','6')
title('Aceleraci\''on $z$ total de cada nodo','interpreter','latex');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Aceleraci\''on total en $z$ [m/s$^2$]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;


%% 2.1 Fijación de nodos para pasar a problema cuasiestático
p = [4 1 0
    4 2 0
    4 3 0
    3 3 0
    3 2 0
    5 3 0];
%%
a
%%
F_vec = zeros(12,3,length(t_vec_full));
for i=1:length(t_vec_full)
    if t_vec_full(i)<=tg
        j=i;
        F_vec(:,:,i) = [1 3 Rz_func(t_vec_full(i))-(g+a_vec(3,i,1))*mPIC/2
                        2 3 Rz_func(t_vec_full(i))-(g+a_vec(3,i,2))*mPIC/2
                        1 1 Rx_func(t_vec_full(i))-a_vec(1,i,1)*mPIC/2
                        2 1 Rx_func(t_vec_full(i))-a_vec(1,i,2)*mPIC/2
                        3 3 L_func(t_vec_full(i))/3-(g+a_vec(3,i,3))*mTela/3
                        4 3 L_func(t_vec_full(i))/3-(g+a_vec(3,i,4))*mTela/3
                        5 3 L_func(t_vec_full(i))/6-(g+a_vec(3,i,5))*mTela/6
                        6 3 L_func(t_vec_full(i))/6-(g+a_vec(3,i,6))*mTela/6
                        3 1 -D_func(t_vec_full(i))/3-a_vec(1,i,3)*mTela/3
                        4 1 -D_func(t_vec_full(i))/3-a_vec(1,i,4)*mTela/3
                        5 1 -D_func(t_vec_full(i))/6-a_vec(1,i,5)*mTela/6
                        6 1 -D_func(t_vec_full(i))/6-a_vec(1,i,6)*mTela/6];
           
   
    elseif t_vec_full(i)<=tr
        F_vec(:,:,i) = [1 3 Rz_vec(i-j)-(g+a_vec(3,i,1))*mPIC/2
                        2 3 Rz_vec(i-j)-(g+a_vec(3,i,2))*mPIC/2
                        1 1 Rx_vec(i-j)-a_vec(1,i,1)*mPIC/2
                        2 1 Rx_vec(i-j)-a_vec(1,i,2)*mPIC/2
                        3 3 L_vec(i-j)/3+Delta_Lg_func(t_vec_full(i))/3-(g+a_vec(3,i,3))*mTela/3
                        4 3 L_vec(i-j)/3+Delta_Lg_func(t_vec_full(i))/3-(g+a_vec(3,i,4))*mTela/3
                        5 3 L_vec(i-j)/6+Delta_Lg_func(t_vec_full(i))/6-(g+a_vec(3,i,5))*mTela/3
                        6 3 L_vec(i-j)/6+Delta_Lg_func(t_vec_full(i))/6-(g+a_vec(3,i,6))*mTela/3
                        3 1 -D_vec(i-j)/3-a_vec(1,i,3)*mTela/3
                        4 1 -D_vec(i-j)/3-a_vec(1,i,4)*mTela/3
                        5 1 -D_vec(i-j)/6-a_vec(1,i,5)*mTela/3
                        6 1 -D_vec(i-j)/6-a_vec(1,i,6)*mTela/3];
 
    elseif t_vec_full(i)<=tg+Delta_tg
        F_vec(:,:,i) = [1 3 Rz_vec(i-j)+Delta_R_func(t_vec_full(i))/2-(g+a_vec(3,i,1))*mPIC/2
                        2 3 Rz_vec(i-j)+Delta_R_func(t_vec_full(i))/2-(g+a_vec(3,i,2))*mPIC/2
                        1 1 Rx_vec(i-j)-a_vec(1,i,1)*mPIC/2
                        2 1 Rx_vec(i-j)-a_vec(1,i,2)*mPIC/2
                        3 3 L_vec(i-j)/3+Delta_Lg_func(t_vec_full(i))/3-(g+a_vec(3,i,3))*mTela/3
                        4 3 L_vec(i-j)/3+Delta_Lg_func(t_vec_full(i))/3-(g+a_vec(3,i,4))*mTela/3
                        5 3 L_vec(i-j)/6+Delta_Lg_func(t_vec_full(i))/6-(g+a_vec(3,i,5))*mTela/6
                        6 3 L_vec(i-j)/6+Delta_Lg_func(t_vec_full(i))/6-(g+a_vec(3,i,6))*mTela/6
                        3 1 -D_vec(i-j)/3-a_vec(1,i,3)*mTela/2
                        4 1 -D_vec(i-j)/3-a_vec(1,i,4)*mTela/2
                        5 1 -D_vec(i-j)/6-a_vec(1,i,5)*mTela/2
                        6 1 -D_vec(i-j)/6-a_vec(1,i,6)*mTela/2];

    elseif t_vec_full(i)<=tr+Delta_tr
        F_vec(:,:,i) = [1 3 Rz_vec(i-j)+Delta_R_func(t_vec_full(i))/2-(g+a_vec(3,i,1))*mPIC/2
                        2 3 Rz_vec(i-j)+Delta_R_func(t_vec_full(i))/2-(g+a_vec(3,i,2))*mPIC/2
                        1 1 Rx_vec(i-j)-a_vec(1,i,1)*mPIC/2
                        2 1 Rx_vec(i-j)-a_vec(1,i,2)*mPIC/2
                        3 3 L_vec(i-j)/3-(g+a_vec(3,i,3))*mTela/3
                        4 3 L_vec(i-j)/3-(g+a_vec(3,i,4))*mTela/3
                        5 3 L_vec(i-j)/6-(g+a_vec(3,i,5))*mTela/6
                        6 3 L_vec(i-j)/6-(g+a_vec(3,i,6))*mTela/6
                        3 1 -D_vec(i-j)/3-a_vec(1,i,3)*mTela/2
                        4 1 -D_vec(i-j)/3-a_vec(1,i,4)*mTela/2
                        5 1 -D_vec(i-j)/6-a_vec(1,i,5)*mTela/2
                        6 1 -D_vec(i-j)/6-a_vec(1,i,6)*mTela/2];

    else
        F_vec(:,:,i) = [1 3 Rz_vec(i-j)-(g+a_vec(3,i,1))*mPIC/2
                        2 3 Rz_vec(i-j)-(g+a_vec(3,i,2))*mPIC/2
                        1 1 Rx_vec(i-j)-a_vec(1,i,1)*mPIC/2
                        2 1 Rx_vec(i-j)-a_vec(1,i,2)*mPIC/2
                        3 3 L_vec(i-j)/3-(g+a_vec(3,i,3))*mTela/3
                        4 3 L_vec(i-j)/3-(g+a_vec(3,i,4))*mTela/3
                        5 3 L_vec(i-j)/6-(g+a_vec(3,i,5))*mTela/6
                        6 3 L_vec(i-j)/6-(g+a_vec(3,i,6))*mTela/6
                        3 1 -D_vec(i-j)/3-a_vec(1,i,3)*mTela/2
                        4 1 -D_vec(i-j)/3-a_vec(1,i,4)*mTela/2
                        5 1 -D_vec(i-j)/6-a_vec(1,i,5)*mTela/2
                        6 1 -D_vec(i-j)/6-a_vec(1,i,6)*mTela/2];
    end
end



%% 2.2 Resolución de la estructura para cada tiempo
u_vec = zeros(data.ni*data.nnod,length(t_vec_full));
r_vec = zeros(data.ni*data.nnod,length(t_vec_full));
sig_vec = zeros(data.nel,length(t_vec_full));
f_vec = zeros(data.nnod*data.ni,length(t_vec_full));
fel_vec = zeros(data.nne*data.ni,data.nel,length(t_vec_full));
% x_new_vec = zeros([size(x),length(t_vec_full)]); % geometria en cada instante de tiempo. x_new_vec[gdl][tiempo]
% x_new_vec(:,:,1) = x;

tic
for i=1:length(t_vec_full)
    a = reshape(a_vec(:,i,:),data.ni,data.nnod)'; % a[nodos,gdl]
    [u_vec(:,i), r_vec(:,i), sig_vec(:,i),f_vec(:,i),fel_vec(:,:,i),l_vec] = solveStructure(data,x,Tn,Td,m,Tm,p,F_vec(:,:,i),a,zeros(data.nel,1));
    % deformacion = reshape(u_vec(:,i),size(x'))';
    % x_new_vec(:,:,i+1) = x_new_vec(:,:,i) + deformacion;
    
end
toc


%%
plot(t_vec_full,r_vec([9,end-4,end],:))
%%
plot3DBars(x,Tn,100,t_vec_full,u_vec,sig_vec,'Pa')
%%

FF_vec = zeros(data.nnod*data.ni,length(t_vec_full));

for i=1:length(t_vec_full)
for j=1:12
    pos=(F_vec(j,1,i)-1)*3+F_vec(j,2,i);
    FF_vec(pos,i) = F_vec(j,3,i);
end
end


%% Safety factor for traction
sf_traction_bar_vec = sigmaYbar*1e6./sig_vec(1:6,:);
sf_traction_cable_vec = sigmaYcable*1e6./sig_vec(7:end,:);

sf_traction_bar_min = min(min(sf_traction_bar_vec(sf_traction_bar_vec>0)));
sf_traction_cable_min = min(min(sf_traction_cable_vec(sf_traction_cable_vec>0)));


%% Safety factor for buckling
sigma_cr_buckle_bar_vec = (pi^2*E1*pi*(2.5e-3/2)^4/4)./(l_vec(1:6).^2*A1)
sigma_cr_buckle_cable_vec = (pi^2*E2*pi*((14e-3/2)^4-(10e-3/2)^4)/4)./(l_vec(7:end).^2*A2)

sf_buckle_bar_vec = sigma_cr_buckle_bar_vec./sig_vec(1:6,:);
sf_buckle_cable_vec = sigma_cr_buckle_cable_vec./sig_vec(7:end,:);

sf_buckle_bar_min = min(min(abs(sf_buckle_bar_vec(sig_vec(1:6,:)<0))));
sf_buckle_cable_min = min(min(abs(sf_buckle_bar_vec(sig_vec(7:end,:)<0))));

%% Gráfico aceleración nodos
for i=1:length(t_vec_full)
    for j=1:6
        a_mod(j,i) = norm(a_vec(:,i,j));
    end
end

figure('Name','Aceleración total de cada nodo');
plot(t_vec_full,a_mod)
xlim([2.5 7])
ylim([-0 .4])
legend('1','2','3','4','5','6')
title("Aceleraci\'on total de cada nodo",'interpreter','latex');
xlabel('Tiempo [s]','interpreter','latex');
ylabel("Aceleraci\'on total [m/s$^2$]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;
