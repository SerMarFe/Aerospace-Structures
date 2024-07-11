function [q, tau, s] = getTangentialStressDistribution(x,Tn,m,Tm,Xc,Ixx,Iyy,Ixy,Sx_ext,Sy_ext)
    % tau is the tangential tension at each node of each element and s is the arc length distance from the
    % first node to the each node of the current element
    nel = length(Tn);
    q = zeros(2, nel);
    tau = zeros(2, nel);
    s = zeros(2, nel);

    % Equivalent inertia and shear:
    Ix = Ixx-Ixy^2/Iyy;
    Iy = Iyy - Ixy^2/Ixx;
    Sx = Sx_ext - Sy_ext*Ixy/Ixx;
    Sy = Sy_ext - Sx_ext*Ixy/Iyy;
    
    % Initial shear force:
    q(1,1) = 0;  % Open section

    for e=1:nel
        t = m(Tm(e), 1);
        dX = x(Tn(e,2),:) - x(Tn(e,1),:); % la majúscula és per indicar que és vector
        l_e = sqrt(sum(dX.^2));
        if e>1
            s(1,e) = s(2,e-1);
            q(1,e) = q(2,e-1);
        end
        s(2,e) = s(1,e) + l_e;
        % Compute shear stress corresponding to first node:
        tau(1, e) = q(1,e)/t;
        % Update shear flow:
        q(2,e) = q(1,e) - Sx*t*l_e*(dX(1)/2 + x(Tn(e,1),1) - Xc(1))/Iy - Sy*t*l_e*(dX(2)/2 + x(Tn(e,1),2) - Xc(2))/Ix;
        % Compute shear stress corresponding to second node:
        tau(2, e) = q(2,e)/t;
    end
end