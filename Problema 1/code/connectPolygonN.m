function tn = connectPolygonN(n)
    tn = zeros(2*n,2);
    for i=1:(n-1)
        tn(i,:) = [i i+1];
        tn(i+n,:) = [n+1 i];
    end
    tn(n,:) = [n 1];
    tn(2*n,:) = [n+1 n];
end