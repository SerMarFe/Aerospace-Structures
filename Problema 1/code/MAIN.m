%% STRUCTURAL PROBLEM CODE STRUCTURE

clear
close all
format long
addpath 'Task 2'/
addpath 'Problema 1/code'

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node

% 1.2 Build geometry (mesh)
% Nodal coordinates matrix
x = [% column_1 = x-coord , column_2 = y-coord , ...    
        0        0
    0.459   -0.054
    1.125        0
    0.315    0.486
    0.864    0.486
];
data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom

% Nodal connectivities matrix
Tn = [% column_1 = element node 1 , column_2 = element node 2, ...
    1   2
    1   4
    2   4
    2   5
    3   5
    4   5
];
data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);

% Material properties matrix
E=71e9;
A1=pi*(37.5^2-34.5^2)*1e-6/4;
A2=pi*(31.2^2-28.8^2)*1e-6/4;
A3=pi*(21^2-19^2)*1e-6/4;
m = [% Each column corresponds to a material property (area, Young's modulus, etc.)
    E    A1    0
    E    A2    0
    E    A3    0
];

% Material connectivities matrix
Tm = [% Each row is the material (row number in 'm') associated to each element
    3
    3
    2
    1
    2
    1
];


% 1.3 Input boundary conditions
% Fixed nodes matrix
p = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    1   1   0
    1   2   0
    3   1   0
    3   2   0
];

% Point loads matrix
W = 75*9.8;
force=75*2.5; % N
F = [% Each row is a point force component | column_1 = node, column_2 = direction (1 = x-direction, 2 = y-direction), column_3 = force magnitude
    2   2   -.45*W
    4   2   -.5*W
    5   1   force
    5   2   -.05*W
];

%% 2) SOLVER

% 2.1.1 Compute element stiffness matrices
Kel = stiffnessFunction(data,x,Tn,m,Tm);

% 2.1.2 Compute element force vectors
fel = forceFunction(data,x,Tn,m,Tm); 

% 2.2 Assemble global stiffness matrix

% Task 2.6
computer = GlobalStiffnessMatrixComputer(Kel,Td,data.ndof); % create an instance of the class
computer.compute() % compute the stiffness matrix using the class method
K = computer.K; % save the stiffness matrix computed using the class
[K1,f] = assemblyFunction(data,Td,Kel,fel); % compute the stiffness matrix using the original code
% Comparison between the stiffness matrix calculated using the class and
% using the original code
if all(all(u1 == u))
    disp("Same result! :)") % Task 2.6
else
    disp("Not same result :(")
end

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data,p);

% 2.3.2 Apply point loads
f = pointLoads(data,f,F); %pointLoads(data,Td,f,F) Td???

%%
% 2.4 Solve system
% Task 2.7

u = zeros(data.ndof,1);
u(vp) = up;
vf = setdiff((1:data.ndof)',vp);

mat = K(vf,vf);
vec = f(vf)-K(vf,vp)*up;
solver = Solver(mat,vec);

[u1,r] = solveSystem(data,K,f,up,vp);

% Comparison between the displacement vector calculated using the class and
% using the original code

u(vf) = solver.solve(0); % solving using the direct solver

if all(all(u1 == u))
    disp("Same result! :)") % Task 2.7
else
    disp("Not same result :(")
end

u(vf) = solver.solve(1); % solving using the iterative solver

if all(all(round(u1,14) == round(u,14)))
    disp("Same result! :)") % Task 2.7
else
    disp("Not same result :(")
end

%%

% 2.5 Compute stress
sig = stressFunction(data,x,Tn,m,Tm,Td,u);

%% 3) POSTPROCESS

scale = 400; % Set a number to visualize deformed structure properly
units = 'Pa'; % Define in which units you're providing the stress vector

plot2DBars(data,x,Tn,u,sig,scale,units);