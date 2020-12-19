function [xkappa] = compute_xkappa(rkappa,rdkappa,l,dl,ddl)

denominator = (dl * dl + (1 - l * rkappa) * (1 - l * rkappa));
if abs(denominator) < 1e-8 
    xkappa = 0;
    return
end
 
denominator = denominator^(1.5);
numerator = rkappa + ddl - 2 * l * rkappa * rkappa - l * ddl * rkappa + l * l * rkappa * rkappa * rkappa + l * dl * rdkappa + 2 * dl * dl * rkappa;
xkappa = numerator / denominator;

end

