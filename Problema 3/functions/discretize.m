function [x_new,Tn] = discretize(x,n)
% Discretiza la geometría, haciendo que cada barra pase a ser n barras
    x_new = zeros((length(x)-1)*n,2);
    Tn = zeros((length(x)-1)*n,2);
    for e = 1:length(x)-1
        for i = 1:n
            dl = (x(e+1,:)-x(e,:))/n;
            x_new(n*(e-1)+i+1,:) = x_new(n*(e-1)+i,:) + dl;
            Tn(n*(e-1)+i,:) = [n*(e-1)+i n*(e-1)+i+1];
        end
    end
end