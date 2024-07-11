function Td = connectDOF(data,Tn)
    nrows = data.nel;
    ncols = data.nne*data.ni;
    Td = zeros(nrows,ncols);
    
    for i=1:nrows
        for j=1:ncols
            col_actual = fix(j/(data.ni+1))+1;
            direccion = mod(j,(data.ni+1))+fix(j/(data.ni+1));
            Td(i,j) = nod2dof(data.ni,Tn(i,col_actual),direccion);
        end
    end
end