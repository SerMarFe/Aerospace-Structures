function [K,f] = assemblyFunction(data,Td,Kel,fel)
    K = zeros(data.ni*data.nnod);
    f = zeros(data.ni*data.nnod,1);
    for e=1:data.nel
        K(Td(e,:),Td(e,:)) = K(Td(e,:),Td(e,:)) + Kel(:,:,e);
        f(Td(e,:),1) = f(Td(e,:),1) + fel(:,e);
    end
end