clear; clc; close all;
addpath('functions\')

% Datos geométricos
c = 2; % [m]
d = .3*c;
h1 = .25*c;
h2 = .15*c;
t1 = 22e-3;
t2 = 15e-3;
t3 = 3.5e-3;
n_subBars = 10;

% Datos del material
E = 210e9;
G = 80e9;

% Geometría y propiedades
x = [0 0
     0 h1/2
     -d h2/2
     -d -h2/2
     0 -h1/2
     0 0];

[x_disc, Tn] = discretize(x,n_subBars);

nel = length(Tn);

m = [t1 E G
     t2 E G
     t3 E G];

Tm = [];
for e = 1:length(x)-1
    if e==1 || e==5
        Tm = [Tm;ones(n_subBars,1)*1];
    end
    if e==2 || e==4
        Tm = [Tm;ones(n_subBars,1)*3];
    end
    if e==3
        Tm = [Tm;ones(n_subBars,1)*2];
    end
end

%% A) Cross-section analysis

% a.1)
[Xc,Xs,Atot,Ixx,Iyy,Ixy,J] = getSectionProperties(x_disc,Tn,m,Tm);
[~,~,~,~,~,~,J_closed,Ain] = getSectionPropertiesClosed(x_disc,Tn,m,Tm);

% a.2) i a.3)
[sigma, s_n] = getNormalStressDistribution(x_disc,Tn,Xc,Ixx,Iyy,Ixy,-1,0);

[q, tau, s] = getTangentialStressDistribution(x_disc,Tn,m,Tm,Xc,Ixx,Iyy,Ixy,0,1);
[q_closed, tau_closed, s_closed] = getTangentialStressDistributionClosed(x_disc,Tn,m,Tm,Xc,Ixx,Iyy,Ixy,0,1,Xs,Ain); % assuming Xs=Xc

[tau_t, s_t] = getTangentialStressDistributionTorsion(x_disc,Tn,m,Tm,J,1);
[tau_t_closed, s_t_closed] = getTangentialStressDistributionTorsionClosed(x_disc,Tn,m,Tm,1,Ain);

%
section_data.x_disc = x_disc;
section_data.Tn = Tn;
section_data.Tm = Tm;
section_data.m = m;
section_data.Xc = Xc;
section_data.Xs = Xs;
section_data.Ixx = Ixx;
section_data.Iyy = Iyy;
section_data.Ixy = Ixy;
section_data.Ain = Ain;
section_data.J = J;
section_data.J_closed = J_closed;

save('section_data.mat','section_data')

%% Gràfics

% a.0) Discretization
hold on;
box on;
grid on; grid minor;
plot(x_disc(:,1),x_disc(:,2),'b-')
plot(x_disc(:,1),x_disc(:,2),'*','color','g')
plot(Xc(1),Xc(2),'*m')
plot(Xs(1),Xs(2),'*c')
legend('Elements','Nodes','Centroid ($\approx$ shear center if closed)','Shear center (opened)','location','northwest','interpreter','latex')
title("Geometry of the wing's section");
xlabel('$-z$ [m]','interpreter','latex');
ylabel("$y$ [m]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1)
xlim([min(x(:,1))-0.6,max(x(:,1))+0.1]);
text(x_disc(1:end-1,1),x_disc(1:end-1,2),sprintfc(' %d',1:numel(x_disc(:,1))-1));
hold off;

%%
% a.1)
% Normal stress (opened section)
figure('Name','Normal stress distribution');
plot(s_n, sigma,'b-')
title('Normal distribution');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("Normal stress [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;


% a.2)
% Shear flow (opened section)
figure('Name','Shear flow distribution (opened)');
plot(s, q,'r-')
title('\textbf{Shear flow distribution (opened)}','Interpreter','latex');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("$q^S_o$ [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

% Tangential stress (opened section)
figure('Name','Tangential stress distribution (opened)');
plot(s, tau,'r-')
title('\textbf{Tangential stress distribution (opened)}','Interpreter','latex');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("$\tau^S_o$ [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

% Torsion tangential stress (opened section)
figure('Name','Torsion tangential stress distribution (opened)');
plot(s, tau_t,'r-')
title('\textbf{Torsion tangential stress distribution (opened)}','Interpreter','latex');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("$\tau^T_o$ [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;



% a.3)
% Shear flow (closed section)
figure('Name','Shear flow distribution (closed)');
plot(s_closed, q_closed,'g-')
title('\textbf{Shear flow distribution (closed)}','Interpreter','latex');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("$q^S_c$ [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

% Tangential stress (closed section)
figure('Name','Tangential stress distribution (closed)');
plot(s_closed, tau_closed,'g-')
title('\textbf{Tangential stress distribution (closed)}','Interpreter','latex');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("$\tau^S_c$ [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

% Torsion tangential stress (closed section)
figure('Name','Tangential stress distribution (closed)');
plot(s_closed, tau_t_closed,'g-')
title('\textbf{Torsion tangential stress distribution (closed)}','Interpreter','latex');
xlabel('Section path ($s$) [m]','interpreter','latex');
ylabel("$\tau^T_c$ [Pa]",'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;