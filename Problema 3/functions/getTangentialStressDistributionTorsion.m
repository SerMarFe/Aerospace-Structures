function [tau_t, s] = getTangentialStressDistributionTorsion(x, Tn, m, Tm, J, Mz_ext)
    % tau is the tangential tension due to torsion at each node of each element and s is the arc length distance from the
    % first node to the each node of the current element
    nel = length(Tn);
    tau_t = zeros(2, nel);
    s = zeros(2, nel);

    % Loop:
    for e= 1:nel
        % Vector element thickness:
        t = m(Tm(e), 1);
        dX = x(Tn(e,2),:) - x(Tn(e,1),:); % la majúscula és per indicar que és vector
        l_e = sqrt(sum(dX.^2));
        if e>1
            s(1,e) = s(2,e-1);
        end
        s(2,e) = s(1,e) + l_e;
        % For open cross-section:
        tau_t(1,e) = Mz_ext*t/J;
        tau_t(2, e) = Mz_ext*t/J;
    end
end