function mNodos = masaNodos(data,Mel,Tn)
    mNodos = zeros(data.nnod,1);
    for e=1:data.nel
        pos = Tn(e,:);
        mNodos(pos,1) = mNodos(pos,1) + Mel(e)/data.nne;
    end
end