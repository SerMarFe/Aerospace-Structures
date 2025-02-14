function sig = stressFunction(data,x,Tn,m,Tm,Td,u)
sig = zeros(data.nel,1);

for e=1:data.nel
    xel = [x(Tn(e,:),:)];
    Eel = m(Tm(e),1);
    sigmael0 = m(Tm(e),3);
    l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2);
    c = (xel(2,1)-xel(1,1))/l;
    s = (xel(2,2)-xel(1,2))/l;
    R = [c s 0 0
         -s c 0 0
         0 0 c s
         0 0 -s c];
    uel = u(Td(e,:));
    sig(e) = (Eel/l)*[-1 0 1 0]*R*uel + sigmael0;
end