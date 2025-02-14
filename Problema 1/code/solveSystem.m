function [u,r] = solveSystem(data,K,f,up,vp)
    u = zeros(data.ndof,1);
    r = zeros(data.ndof,1);

    u(vp) = up;
    vf = setdiff((1:data.ndof)',vp);
    
    u(vf) = K(vf,vf)\(f(vf)-K(vf,vp)*up);
    r(vp) = K(vp,:)*(u) - f(vp);

end