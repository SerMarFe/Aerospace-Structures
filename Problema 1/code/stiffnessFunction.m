function Kel = stiffnessFunction(data,x,Tn,m,Tm)

    Kel = zeros(4,4,data.nel);
 
    for e=1:data.nel
    xel = [x(Tn(e,:),:)];
    Eel = m(Tm(e),1);
    Ael = m(Tm(e),2);
    
    l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2);
    c = (xel(2,1)-xel(1,1))/l;
    s = (xel(2,2)-xel(1,2))/l;
    
    Kel(:,:,e) = (Eel*Ael/l)*[
        c^2 c*s -c^2 -c*s
        c*s s^2 -c*s -s^2
        -c^2 -c*s c^2 c*s
        -c*s -s^2 c*s s^2
        ];
    end

end