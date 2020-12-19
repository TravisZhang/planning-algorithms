syms dx dy d2x d2y d3x d3y
a = dx * d2y - dy * d2x;
b = dx * d3y - dy * d3x;
c = dx * d2x + dy * d2y;
d = dx^2 + dy^2;
dkappa = (b * d - 3 * a * c) / d^3;
dkappa = expand(dkappa)