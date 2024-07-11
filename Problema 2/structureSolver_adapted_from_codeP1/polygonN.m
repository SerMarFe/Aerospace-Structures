function x = polygonN(n,R)
    x = zeros(n,2);
    a = (360/n);
    for i=1:n
        cordx = R*cosd(a*(i-1) + a/2 - (90 + a));
        cordy = R*sind(a*(i-1) + a/2 - (90 + a));
        x(i,:) = [cordx cordy];
    end
end