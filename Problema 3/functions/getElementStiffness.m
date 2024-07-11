function [Kel] = getElementStiffness(e, x, Tn, m, Tm)
% Function in order to obtain the element stiffness matrix
    % Vector element thickness:
    E_e = m(Tm(e), 1);
    G_e = m(Tm(e), 2);
    I_e = m(Tm(e), 3);
    J_e = m(Tm(e), 4);
    l_e = abs(x(Tn(e,2)) - x(Tn(e,1))); % la majúscula és per indicar que és vector
    % Compute element stiffness matrix for bending:
    Kb = (E_e*I_e/l_e^3)*[
                           12 6*l_e 0 -12 6*l_e 0;
                           6*l_e 4*l_e^2 0 -6*l_e 2*l_e^2 0;
                           0 0 0 0 0 0;
                           -12 -6*l_e 0 12 -6*l_e 0;
                           6*l_e 2*l_e^2 0 -6*l_e 4*l_e^2 0;
                           0 0 0 0 0 0];
    % Compute element stiffness matrix for torsion:
    Kt = (G_e*J_e/l_e)*[
                        0 0 0 0 0 0;
                        0 0 0 0 0 0;
                        0 0 1 0 0 -1;
                        0 0 0 0 0 0;
                        0 0 0 0 0 0;
                        0 0 -1 0 0 1];
    % Compute total element stiffness matrix:
    Kel = Kb + Kt;
end