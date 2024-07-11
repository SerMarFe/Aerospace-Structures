clear; clc; close all;

addpath('functions\')
load('section_data.mat')

g = 9.81; % [m/s^2]
% Datos geométricos
chi_s = -section_data.Xc(1); % [m] shear center caso cerrado (approx= centroide)
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
dx = b/512;
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
J = section_data.J_closed; % [m^4] Módulo de torsión (polar inertia) (caso abierto)
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


%% B) Beam analysis (closed)
[u,r, x_el,S_el,Mb_el,Mt_el] = solveStructure(data,Tn,Td,m,Tm,p,distributed_load,distributed_torsion,F);

u_3rows = reshape(u,3,[]);


%% Gràfics

% b.1)
figure('Name','Closed section');
t = tiledlayout(3,1);
title(t,'\textbf{Closed section}','interpreter','latex')

% Vertical deflection (closed section)
ax1 = nexttile;
plot(ax1,x,u_3rows(1,:)*100,'b-')
title(ax1,'\textbf{Vertical deflection}','interpreter','latex');
ylabel("$y$ [cm]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;
%daspect([1 4 .5])

% Bending rotation (closed section)
ax2 = nexttile;
plot(ax2,x,u_3rows(2,:)*180/pi,'g-')
title(ax2,'\textbf{Bending rotation}','interpreter','latex');
ylabel("$\theta_z$ [deg]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

% Torsion rotation (closed section)
ax3 = nexttile;
plot(ax3,x,u_3rows(3,:)*180/pi,'r-')
title(ax3,'\textbf{Torsion rotation}','interpreter','latex');
ylabel("$\theta_x$ [deg]",'Interpreter','latex');
xlabel('Wingspan position [m]','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

%

figure('Name','Closed section');
t = tiledlayout(3,1);
title(t,'\textbf{Closed section}','interpreter','latex')

% Shear force (closed section)
ax1 = nexttile;
plot(ax1,x_el,S_el/1000,'b-')
title(ax1,'\textbf{Shear force}','interpreter','latex');
ylabel("$S_y$ [kN]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;
%daspect([1 4 .5])

% Bending moment (closed section)
ax2 = nexttile;
plot(ax2,x_el,Mb_el/1000,'g-')
title(ax2,'\textbf{Bending moment}','interpreter','latex');
ylabel("$M_z$ [kN\,m]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

% Torsion moment (closed section)
ax3 = nexttile;
plot(ax3,x_el,Mt_el/1000,'r-')
title(ax3,'\textbf{Torsion moment}','interpreter','latex');
ylabel("$M_x$ [kN\,m]",'Interpreter','latex');
xlabel('Wingspan position [m]','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;


%% B2 Study of convergence

[discretization_factor, tip_deflection_error, tip_torsion_error] = getTipDeflectoinAndTorsionClosedError(512);

% Gráfico de la convergencia
hold on;
box on;
grid on; grid minor;
plot(discretization_factor,tip_deflection_error*100,'bx:','LineWidth',1)
plot(discretization_factor,tip_torsion_error*100,'mx:','LineWidth',1)
set(gca,'xscale','log')

ticks = [1, 2, 4, 8, 16, 32, 64, 128];
tickLabels = cellstr(num2str(ticks'));
set(gca, 'XTick', ticks, 'XTickLabel', tickLabels);

legend('Deflection','Twist','location','best','interpreter','latex')
title("\textbf{Study of Convergence (closed)}",'interpreter','latex');
xlabel('Discretization factor ($n$) [$h_{el} = b/n$]','interpreter','latex');
ylabel("Relative error (\%)",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1)
hold off;



%% C) Most critical point (Von Misses)

aplana = @(vec)[vec(1,1) vec(1,2:end) vec(2,end)];

cross_section_max_sigma_vm_vec = zeros(1,data.nel);
node_cross_section_max_sigma_vm_vec = zeros(1,data.nel);
for i=1:data.nel
    [~,tau,~] = getTangentialStressDistributionClosed(section_data.x_disc,section_data.Tn,section_data.m,section_data.Tm,section_data.Xc,section_data.Ixx,section_data.Iyy,section_data.Ixy,0,S_el(1,i),section_data.Xc,section_data.Ain);
    [tau_t,~] = getTangentialStressDistributionTorsionClosed(section_data.x_disc,section_data.Tn,section_data.m,section_data.Tm,-Mt_el(1,i),section_data.Ain);
    [sigma,~] = getNormalStressDistribution(section_data.x_disc,section_data.Tn,section_data.Xc,section_data.Ixx,section_data.Iyy,section_data.Ixy,Mb_el(1,i),0);
    tau = aplana(tau); % tau de cada nodo (cada posición es el nodo correspondiente de la discretización de la sección transversal)
    tau_t = aplana(tau_t); % tau_t de cada nodo (cada posición es el nodo correspondiente de la discretización de la sección transversal)
    sigma = aplana(sigma); % sigma de cada nodo (cada posición es el nodo correspondiente de la discretización de la sección transversal)

    sigmas_vm = sqrt(sigma.^2+3*(tau+tau_t).^2);
    [cross_section_max_sigma_vm_vec(i), node_cross_section_max_sigma_vm_vec(i)] = max(sigmas_vm);
end

[wingspan_max_sigma_vm, node_wingspan_max_sigma_vm] = max(cross_section_max_sigma_vm_vec)
node_cross_section_max_sigma_vm = node_cross_section_max_sigma_vm_vec(node_wingspan_max_sigma_vm)

position_wing_max_sigma_vm = (node_wingspan_max_sigma_vm-1)*dx
position_cross_section_max_sigma_vm = (node_cross_section_max_sigma_vm-1)*norm(diff(section_data.x_disc(1:2,:)))

% Cross section maximum sigma Von Misses
hold on;
box on;
grid on; grid minor;
plot(section_data.x_disc(:,1),section_data.x_disc(:,2),'b-')
plot(section_data.x_disc(:,1),section_data.x_disc(:,2),'*','color','g')
plot(section_data.x_disc(node_cross_section_max_sigma_vm,1),section_data.x_disc(node_cross_section_max_sigma_vm,2),'*','color','r')
legend('Elements','Nodes','Max $\sigma_{vm}$','location','best','interpreter','latex')
title(sprintf("Most demanded node (closed, section at %.2f m)",position_wing_max_sigma_vm));
xlabel('$-z$ [m]','interpreter','latex');
ylabel("$y$ [m]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1)
%xlim([min(x(:,1))-0.1,max(x(:,1))+0.1]);
text(section_data.x_disc(1:end-1,1),section_data.x_disc(1:end-1,2),sprintfc(' %d',1:numel(section_data.x_disc(:,1))-1));
hold off;