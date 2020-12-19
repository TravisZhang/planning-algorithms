syms a0 a1 a2
syms b0 b1 b2
syms c00 c01 c02 c10 c11 c12 c20 c21 c22
a = [a0 a1 a2].';
b = [b0 b1 b2].';
C = [c00 c01 c02; c10 c11 c12; c20 c21 c22]; % C not symmetric
C = [c00 c01 c02; c01 c11 c12; c02 c12 c22]; % C symmetric

D = a.'*C*b;
E = b.'*C*a;
D = expand(D)
E = expand(E)
fprintf(1,'D == E: %d\n',D==E);
% Conclusion
% if C not symmetric, D,E are not the same
% if C symmetric, D E are the same !!!