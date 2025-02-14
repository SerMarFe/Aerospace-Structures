function Td = connectDOF(data,Tn)
    nrows = data.nel;
    ncols = data.nne*data.ni;
    Td = zeros(nrows,ncols);
    
    for i=1:nrows
        for j=1:ncols
            Td(i,j) = data.ni*(Tn(i,mod(j,data.nne)+fix(j/data.nne))-1) +mod(j,data.ni+1)+fix(j/(data.ni+1));
        end
    end
end