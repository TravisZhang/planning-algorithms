function [xdkappa] = compute_xdkappa(rkappa,rdkappa,rddkappa,l,dl,ddl,dddl)

kappa = 0;
denominator = (dl * dl + (1 - l * rkappa) * (1 - l * rkappa));
if abs(denominator) > 1e-8 
    denominator = denominator^(1.5);
    numerator = rkappa + ddl - 2 * l * rkappa * rkappa - l * ddl * rkappa + l * l * rkappa * rkappa * rkappa + l * dl * rdkappa + 2 * dl * dl * rkappa;
    kappa = numerator / denominator;
end
 
one_minus_l_rk = 1 - rkappa * l;
one_minus_l_rk_prime = -(rdkappa * l + rkappa * dl);
one_minus_l_rk_dprime = -(ddl * rkappa + 2.0 * dl * rdkappa + l * rddkappa);

tan_theta_diff = dl / one_minus_l_rk;
theta_diff = atan2(dl, one_minus_l_rk);
cos_theta_diff = cos(theta_diff);
cos_theta_diff_2 = cos_theta_diff * cos_theta_diff;
cos_theta_diff_3 = cos_theta_diff_2 * cos_theta_diff;
theta_diff_prime = kappa * one_minus_l_rk / cos_theta_diff - rkappa;

theta_diff_dprime = (dddl - one_minus_l_rk_dprime * tan_theta_diff - 2.0 * one_minus_l_rk_prime * theta_diff_prime / cos_theta_diff_2) * cos_theta_diff_3 / one_minus_l_rk;
theta_diff_dprime = (theta_diff_dprime + 2.0 * theta_diff_prime * theta_diff_prime) / cos_theta_diff;
res = theta_diff_dprime + rdkappa - kappa * (one_minus_l_rk_prime + one_minus_l_rk * tan_theta_diff * theta_diff_prime) / cos_theta_diff;
res = res * cos_theta_diff_2 / (one_minus_l_rk * one_minus_l_rk);

xdkappa = res;

end

