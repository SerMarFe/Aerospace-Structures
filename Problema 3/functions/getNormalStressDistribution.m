function [sigma, s] = getNormalStressDistribution(x,Tn,Xc,Ixx,Iyy,Ixy,Mx_ext,My_ext)
    % sigma is the normal tension at each node of each element and s is the arc length distance from the
    % first node to the each node of the current element
    nel = length(Tn);
    sigma = zeros(2, nel);
    s = zeros(2, nel);
    Ix = Ixx-Ixy^2/Iyy;
    Iy = Iyy - Ixy^2/Ixx;
    Mx = Mx_ext + My_ext*Ixy/Iyy;
    My = My_ext + My_ext*Ixy/Ixx;

    for e=1:nel
        dX = x(Tn(e,2),:) - x(Tn(e,1),:); % la majúscula és per indicar que és vector
        l_e = sqrt(sum(dX.^2));
        if e>1
            s(1,e) = s(2,e-1);
        end
        s(2,e) = s(1,e) + l_e;
        sigma(1,e) = (x(Tn(e,1),2)-Xc(2))*Mx/Ix - (x(Tn(e,1),1)-Xc(1))*My/Iy;
        sigma(2,e) = (x(Tn(e,2),2)-Xc(2))*Mx/Ix - (x(Tn(e,2),1)-Xc(1))*My/Iy;
    end
end