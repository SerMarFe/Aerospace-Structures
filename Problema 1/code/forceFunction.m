function fel = forceFunction(data,x,Tn,m,Tm)

    fel = zeros(data.nne*data.ni,data.nel);
    for e=1:data.nel
        xel = [x(Tn(e,:),:)];
        Ael = m(Tm(e),2);
        sigma0el = m(Tm(e),3);
        
        l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2);
        c = (xel(2,1)-xel(1,1))/l;
        s = (xel(2,2)-xel(1,2))/l;
        R = [
            c s 0 0
            -s c 0 0
            0 0 c s
            0 0 -s c];
        
        fel(:,e) = -sigma0el*Ael*R.'*[-1 0 1 0]';
    end
end