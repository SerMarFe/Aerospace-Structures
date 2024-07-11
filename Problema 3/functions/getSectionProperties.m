function [Xc,Xs,Atot,Ixx,Iyy,Ixy,J] = getSectionProperties(x,Tn,m,Tm)
    % Centroid calculation    
    Xc = [0, 0];
    Atot = 0;
    nel = length(Tn);
    l_e = zeros(nel,1);
    A_e = zeros(nel,1);
    X_e = zeros(nel,2);
    
    for e = 1:nel
        t_e = m(Tm(e),1);
        dX = x(Tn(e,2),:) - x(Tn(e,1),:); % la majúscula és per indicar que és vector
        l_e(e) = sqrt(sum(dX.^2));
        A_e(e) = t_e*l_e(e);
        Atot = Atot + A_e(e);
        
        X_e(e,:) = (x(Tn(e,1),:)+x(Tn(e,2),:))/2;
        Xc = Xc + X_e(e,:)*A_e(e);   
    end
    Xc = Xc/Atot;
    % Xc = A_e'*X_e/Atot; % equivalent
    
    % Inertia calculation
    Ixx = 0;
    Iyy = 0;
    Ixy = 0;
    J = 0;

    for e = 1:nel
        t_e = m(Tm(e),1);
        dX = x(Tn(e,2),:) - x(Tn(e,1),:);
        
        Ixx = Ixx + (1/12)*A_e(e)*dX(2)^2 + A_e(e)*(X_e(e,2)-Xc(2))^2;
        Iyy = Iyy + (1/12)*A_e(e)*dX(1)^2 + A_e(e)*(X_e(e,1)-Xc(1))^2;
        Ixy = Ixy + (1/12)*A_e(e)*dX(1)*dX(2) + A_e(e)*(X_e(e,2)-Xc(2))*(X_e(e,1)-Xc(1));
        J = J + l_e(e)*t_e^3/3;
    end

    % Shear calculation
    Q = [0, 0]; % per quan actua tallant només en x i només en y (per poder trobar la component y i la x del centre de tallant, respectivament) 
    Xs = Xc;

    for e = 1:nel
        t_e = m(Tm(e),1);
        dX = x(Tn(e,2),:) - x(Tn(e,1),:);
        A = [(Iyy*dX(2)/2-Ixy*dX(1)/2)/(Ixx*Iyy-Ixy^2), (Ixx*dX(1)/2-Ixy*dX(2)/2)/(Ixx*Iyy-Ixy^2)];
        B = [(Iyy*(x(Tn(e,1),2)-Xc(2)) - Ixy*(x(Tn(e,1),1)-Xc(1)))/(Ixx*Iyy-Ixy^2), (Ixx*(x(Tn(e,1),1)-Xc(1)) - Ixy*(x(Tn(e,1),2)-Xc(2)))/(Ixx*Iyy-Ixy^2)];
        C = (x(Tn(e,1),1)-Xc(1))*dX(2) - (x(Tn(e,1),2)-Xc(2))*dX(1);
        Xs = Xs + C*(Q + ((A/3+B/2)*t_e*l_e(e).*[-1, 1]));
        Q = Q + t_e*l_e(e)*(A+B).*[-1, 1];
    end
end