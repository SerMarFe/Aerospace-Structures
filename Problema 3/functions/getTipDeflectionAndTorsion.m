function [deflection_tip,torsion_tip] = getTipDeflectionAndTorsion(n)
% 'n' being the factor that discretizes the total length of the beam (wing)
% according to 0:b/n:b

g = 9.81; % [m/s^2]
% Datos geométricos
load('section_data.mat')
chi_s = -section_data.Xs(1); % [m] shear center caso abierto

b = 16; % [m]
be = 0.25*b; % [m]
c = 2; % [m]
ze = 0.3*c; % [m]
za = 0.25*c; % [m]
zm = 0.48*c; % [m]
chi_p = 0.3*c; % [m]
Me = 2100; % [kg]
weight_distribution = 140*g; % [N/m]

vinf = 750/3.6; % [m/s]
rho = 1.225; % [kg/m^3]
Cl = 0.1;
lift_func = @(x)0.5*rho*vinf^2*c*Cl*sqrt(1-(x/b)^2);

% Geometría y propiedades
dx = b/n;
x = 0:dx:b;
data.x = x;
data.nnod = numel(x);
data.nel = numel(x)-1;
data.nne = 2;
data.ni = 3;
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom
Tn = zeros(data.nel,2);
for e=1:data.nel
    Tn(e,:) = [e, e+1];
end
Td = connectDOF(data,Tn);

p = [1 1 0 % no desplazamiento vertical en el empotramiento
     1 2 0 % no giro en el empotramiento
     1 3 0]; % no giro de torsion en el empotramiento

% Datos del material
E = 210e9;
G = 80e9;
I = section_data.Ixx; % [m^4] Inercia en el eje x local de la sección transversal (el eje z global del dibujo)
J = section_data.J; % [m^4] Módulo de torsión (polar inertia) (caso abierto)
m = [E G I J];
Tm = ones(data.nel,1); % Toda el ala presenta la misma sección con mismo material


% Fuerza distribuida equivalente resultante en cada ELEMENTO
distributed_lift_nodes = arrayfun(lift_func,x);
distributed_lift_elements = (distributed_lift_nodes(1:end-1) + distributed_lift_nodes(2:end))/2;
distributed_load = distributed_lift_elements - weight_distribution;

% Momento torsor distribuido equivalente resultante en cada ELEMENTO
distributed_torsion_lift = (chi_s+chi_p-za)*distributed_lift_elements;
distributed_torsion = distributed_torsion_lift - (chi_s+chi_p-zm)*weight_distribution;

% Fuerzas puntuales y momentos/torsores en cada NODO
motor_node_position = find(x==be);
F = [motor_node_position 1 -Me*g
     motor_node_position 3 -Me*g*(chi_s+chi_p-ze)]; % Point loads (only the engine)


%% B) Beam analysis (opened)
[u,~, ~,~,~,~] = solveStructure(data,Tn,Td,m,Tm,p,distributed_load,distributed_torsion,F);
u_3rows = reshape(u,3,[]);
deflection_tip = u_3rows(1,end);
torsion_tip = u_3rows(3,end);