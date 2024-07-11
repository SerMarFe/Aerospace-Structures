function [x_el, S_el, Mb_el, Mt_el] = getElementInternalForces(data, Tn, Td, Kel, u)
% This function allows to compute the loads diagram along the beam. The
% outputs can then be used to get the corresponding stresses distributions
% on the cross-section.

x = data.x;

% Initialize vectors:
x_el = zeros(2, data.nel);
S_el = zeros(2, data.nel);
Mb_el = zeros(2, data.nel);
Mt_el = zeros(2, data.nel);


% Loop:
for e=1:data.nel
    % 1. Retrieve element's vertices coordinates:
    x_el(:,e) = x(Tn(e,:))';

    % 2. Retrieve the DOFs displacements of each node:
    u_el = u(Td(e,:)); % size is (data.nne*data.ni,1)

    % 3. Compute internal forces
    fint_el = Kel(:,:,e)*u_el;

    % 4. Assign cross-section loads
    S_el(1,e) = -fint_el(1);
    S_el(2,e) = fint_el(4);
    Mb_el(1,e) = -fint_el(2);
    Mb_el(2,e) = fint_el(5);
    Mt_el(1,e) = -fint_el(3);
    Mt_el(2,e) = fint_el(6);
end