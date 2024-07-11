function f = pointLoads(data,f,F)
    node_to_dof = @(ni,node,direction) ni*(node-1) + direction;
    Fext = zeros(size(f));
    Fext(node_to_dof(data.ni,F(:,1),F(:,2))) = F(:,3);
    f = f + Fext;
end
