clear; close all;

E2=210e9; % [Pa]
A2=3.8e-6; % [m^2]
I2 = 1.15e-12; % [m^4]
R = 0.35;
sigmaCrEuler = pi^2*E2*I2/(R^2*A2); % [Pa]
SF = 2.5;
sigmaAdm = sigmaCrEuler/SF; % [Pa]


% Configución de las opciones de optimización
tolerance = 1e-3;
options = optimset('Display','on','TolFun', tolerance);

% Punto inicial
initial_guess = 1e7;

func = @(x) abs(2.5 - sigmaCrEuler/getSigmaCompressionMax(x));

% Encuentra la raíz con fminunc o fminsearch
sigma0Good = fminsearch(func, initial_guess, options)
SFactual = sigmaCrEuler/getSigmaCompressionMax(sigma0Good)

%% Per visulitzar la funció a trobar el mínim o l'arrel
%set(0, 'DefaultTextInterpreter', 'latex');

xx=linspace(sigma0Good-1e6,sigma0Good+1e6,4000);
plot(xx,arrayfun(func,xx),'-r','LineWidth',2)
xlabel('Pretensi\''on inicial ($\sigma_0)$ [Pa]','Interpreter','latex')
ylabel('Funci\''on a encontrar la ra\''iz','Interpreter','latex')
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$|2.5 - \frac{\sigma_{u}}{\sigma_{c_{max}}}|$','interpreter','latex')
grid on; grid minor;