function sig = stressFunction(data,x,Tn,m,Tm,Td,u)
sig = zeros(data.nel,1);

for e=1:data.nel
    xel = [x(Tn(e,:),:)];
    Eel = m(Tm(e),1);
    l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2+(xel(2,3)-xel(1,3))^2);
    dx = (xel(2,1)-xel(1,1))/l;
    dy = (xel(2,2)-xel(1,2))/l;
    dz = (xel(2,3)-xel(1,3))/l;
    R = [
        dx dy dz 0 0 0;
        0 0 0 dx dy dz
        ];
    uel = u(Td(e,:));
    sig(e) = (Eel/l)*[-1 1]*R*uel;
end