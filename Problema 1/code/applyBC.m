function [up,vp] = applyBC(data,p)
    node_to_dof = @(ni,node,direction) ni*(node-1) + direction;
    up = p(:,3);
    vp = node_to_dof(data.ni,p(:,1),p(:,2));
end