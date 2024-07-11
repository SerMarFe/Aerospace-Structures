function [u,r, x_el,S_el,Mb_el,Mt_el] = solveStructure(data,Tn,Td,m,Tm,p,distributed_load,distributed_torsion,F)
    x = data.x;

    Kel = zeros(data.ni*data.nne,data.ni*data.nne,data.nel);
    fel = zeros(data.ni*data.nne,data.nel);

    % 1 Compute the element stiffness matrices and element force vectors
    for e=1:data.nel
        Kel(:,:,e) = getElementStiffness(e,x,Tn,m,Tm); % stiffness matrix of each element [desplazamiento en y, giro (theta) en z, giro (por torsion) en x]
        fel(:,e) = getElementForce(e,x,Tn,distributed_load(e),distributed_torsion(e)); % force on the nodes of each element due to the distributed load and distributed torsor moment
    end
    
    % 2 Assemble global stiffness matrix
    [K,f] = assemblyFunction(data,Td,Kel,fel); % spy(K) (ineficient, millor fer servir sparse, etc...)
    
    % 3 Apply point loads
    f = pointLoads(data,f,F);
    
    % 3 Apply prescribed DOFs
    [up,vp] = applyBC(data,p);
    
    % 4 Solve system
    [u,r] = solveSystem(data,K,f,up,vp);
    
    % 5 Get the loads along the beam
    [x_el, S_el, Mb_el, Mt_el] = getElementInternalForces(data,Tn,Td,Kel,u);
end