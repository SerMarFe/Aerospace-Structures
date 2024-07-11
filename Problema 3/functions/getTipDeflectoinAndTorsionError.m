function [n,tip_deflection_error,tip_torsion_error] = getTipDeflectoinAndTorsionError(n_ref)
    [tip_deflection_ref, tip_torsion_ref] = getTipDeflectionAndTorsion(n_ref);
    
    error = @(x,ref)abs((x-ref)/ref);
    
    for i=1:5
        n(i)=2^(i+1);
        [tip_deflection, tip_torsion] = getTipDeflectionAndTorsion(n(i));
        tip_deflection_error(i) = error(tip_deflection,tip_deflection_ref);
        tip_torsion_error(i) = error(tip_torsion,tip_torsion_ref);
    end
end