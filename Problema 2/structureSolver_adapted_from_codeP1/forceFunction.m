function [fel, l_vec] = forceFunction(data,x,Tn,m,Tm,a,sigmas0)
% 'a' contiene la aceleración de cada nodo x,y,z para un instante dado. a[nodo,gdl]
% de todos los instantes: tip: reshape(a_vec(:,temps,:),data.ni,data.nnod))'
    fel = zeros(data.nne*data.ni,data.nel);
    l_vec = zeros(data.nel,1);
    for e=1:data.nel
        xel = [x(Tn(e,:),:)];
        ael = [a(Tn(e,:),:)];
        Ael = m(Tm(e),2);
        rhoel = m(Tm(e),4);
        sigma0el = sigmas0(e);
        % sigma0el = m(Tm(e),3);
        
        l = sqrt((xel(2,1)-xel(1,1))^2+(xel(2,2)-xel(1,2))^2+(xel(2,3)-xel(1,3))^2);
        l_vec(e) = l;
        dx = (xel(2,1)-xel(1,1))/l;
        dy = (xel(2,2)-xel(1,2))/l;
        dz = (xel(2,3)-xel(1,3))/l;
        
        R = [
            dx dy dz 0 0 0;
            0 0 0 dx dy dz
            ];
        
        g = [0;0;-9.81];
        fel(:,e) = (l*rhoel*Ael/2)*[(g-ael(1,:)');
                                  (g-ael(2,:)')] + -sigma0el*Ael*R.'*[-1 1 ]';
    end
end