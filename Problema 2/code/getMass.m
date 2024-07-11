function Mel = getMass(data,x,Tn,m,Tm)
    Mel = zeros(data.nel,1);
    for e=1:data.nel
        xel = [x(Tn(e,:),:)];
        Densityel = m(Tm(e),4);
        Ael = m(Tm(e),2);
        l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2+(xel(2,3)-xel(1,3))^2);   
        Mel(e) = Ael*l*Densityel;
    end
end