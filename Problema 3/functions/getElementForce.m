function [fel] = getElementForce(e, x, Tn, fe, me)
% Element force function, this algorithm for this function will take into account that a uniformly distributed vertical force.
    l_e = abs(x(Tn(e,2)) - x(Tn(e,1))); % la majúscula és per indicar que és vector
    % Contribution of distributed shear load:
    fb = fe*l_e*[1/2; l_e/12; 0; 1/2; -l_e/12; 0];
    % Contribution of distributed torsion:
    ft = me*l_e*[0; 0; 1/2; 0; 0; 1/2];
    % Compute total element force vector:
    fel = fb + ft;
end