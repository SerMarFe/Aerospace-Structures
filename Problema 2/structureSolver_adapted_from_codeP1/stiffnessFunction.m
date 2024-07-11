function Kel = stiffnessFunction(data,x,Tn,m,Tm)

    Kel = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
 
    for e=1:data.nel
    xel = x(Tn(e,:),:);
    Eel = m(Tm(e),1);
    Ael = m(Tm(e),2);
    
    l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2+(xel(2,3)-xel(1,3))^2);
    dx = (xel(2,1)-xel(1,1))/l;
    dy = (xel(2,2)-xel(1,2))/l;
    dz = (xel(2,3)-xel(1,3))/l;
    
    R = [
        dx dy dz 0 0 0;
        0 0 0 dx dy dz
        ];

    K_prima = (Eel*Ael/l)*[
        1 -1
        -1 1
        ];

    Kel(:,:,e) = R.'*K_prima*R;
    end

end