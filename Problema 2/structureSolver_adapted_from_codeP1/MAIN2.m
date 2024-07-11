%% STRUCTURAL PROBLEM CODE STRUCTURE

clear
close all
format long

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node

% 1.2 Build geometry (mesh)
% Nodal coordinates matrix
% column_1 = x-coord , column_2 = y-coord , ...    
n = 100; % rodes com polígon regular de n costats
R = 35e-2; % [m]
d = 1.125; % [m]
x1 = polygonN(n,R);
x = [x1
    0 0];
%    x1(:,1)+d x1(:,2)
%    d 0];

data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom

% Nodal connectivities matrix

% column_1 = element node 1 , column_2 = element node 2, ...
tn = connectPolygonN(n);
Tn = tn;
   %  tn+(n+1)];

data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);

% Material properties matrix
E1=70e9; % [Pa]
E2=210e9; % [Pa]
A1=140e-6; % [m^2]
A2=3.8e-6; % [m^2]
I1 = 1420e-12; % [m^4]
I2 = 1.15e-12; % [m^4]
sigmaCrEuler = pi^2*E2*I2/(R^2*A2); % [Pa]
%sigma0=4.525452414038591e+07; %trasera pretensión
%sigma0 = 1.105280041599181e+08; %delantera pretensión
sigma0=0;
SF = 2.5;
sigmaAdm = sigmaCrEuler/SF; % [Pa]
m = [% Each column corresponds to a material property (area, Young's modulus, etc.) and inertia [m^4]
    E1    A1    0
    E2    A2    sigma0
];


% Material connectivities matrix
tm = ones(n,1);
Tm = [% Each row is the material (row number in 'm') associated to each element
 tm*1
 tm*2
% tm*1
% tm*2
];

% 1.3 Input boundary conditions
% Fixed nodes matrix
p = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    1   1   0
    1   2   0
    2   2   0
   % (n+1)+1   1   0
   % (n+1)+1   2   0
   % (n+1)+2   2   0
];

% Point loads matrix
Rr = 1e2*[-0.011105555555555 3.879300000000002]; % [N]
Rf = 1e2*[-1.863894444444445 3.470700000000001]; % [N]
F = [% Each row is a point force component | column_1 = node, column_2 = direction (1 = x-direction, 2 = y-direction), column_3 = force magnitude
    n+1 1   -1*Rr(1)
    n+1 2   -1*Rr(2)
   % 2*(n+1) 1   -1*Rf(1) 
   % 2*(n+1) 2   -1*Rf(2)
];

%% 2) SOLVER

% 2.1.1 Compute element stiffness matrices
Kel = stiffnessFunction(data,x,Tn,m,Tm);

% 2.1.2 Compute element force vectors
fel = forceFunction(data,x,Tn,m,Tm); 

% 2.2 Assemble global stiffness matrix
[K,f] = assemblyFunction(data,Td,Kel,fel);

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data,p);

% 2.3.2 Apply point loads
f = pointLoads(data,f,F); %pointLoads(data,Td,f,F) Td???

% 2.4 Solve system
[u,r] = solveSystem(data,K,f,up,vp);

% 2.5 Compute stress
sig = stressFunction(data,x,Tn,m,Tm,Td,u);

sigRad = sig((n+1):2*n);
%               (3*n+1):4*n]);

%% 3) POSTPROCESS

scale = 200; % Set a number to visualize deformed structure properly
units = 'Pa'; % Define in which units you're providing the stress vector

plot2DBars(data,x,Tn,u,sig,scale,units);
ylim([min(x(:,2))-0.1,max(x(:,2))+0.1]);

%%
% close all
% hold on;
% title('Numeración de nodos')
% %plot(x(:,1),x(:,2),'red');c
% scatter(x(:,1),x(:,2),'red','filled')
% % make data labels:
% xlim([min(x(:,1))-0.1,max(x(:,1))+0.1]);
% text(x(:,1),x(:,2),sprintfc(' %d',1:numel(x(:,1))))
% axis equal;