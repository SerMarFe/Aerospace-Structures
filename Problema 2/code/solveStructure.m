function [u,r,sig] = solveStructure(data,x,Tn,Td,m,Tm,p,F)
% Resuelve la estructura y devuleve los desplazamientos, reacciones y
% tensiones

%% 2) SOLVER

% 2.1.1 Compute element stiffness matrices
Kel = stiffnessFunction(data,x,Tn,m,Tm);

% 2.1.2 Compute element force vectors
fel = forceFunction(data,x,Tn,m,Tm,ael); 

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

end