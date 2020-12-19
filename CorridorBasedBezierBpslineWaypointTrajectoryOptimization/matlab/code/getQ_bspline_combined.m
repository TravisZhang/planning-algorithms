function [Q] = getQ_bspline_combined(n, p)
% Calculate smoothness cost Js by elastic band method
% n+1 is number of control points P
% p is degree of bspline, p >= 2(must)
% note: because the forces are vectors, we need to combine xy to optimize
% together

% method 1 elastic band method
Q = zeros(2*(n+1),2*(n+1));
Q_temp = [1 -2 1;-2 4 -2;1 -2 1];
for i = p:n-p+2
    Q(i-1:i+1,i-1:i+1) = Q(i-1:i+1,i-1:i+1) + Q_temp;
end
idx_start = n+1;
for i = p:n-p+2
    Q(idx_start+i-1:idx_start+i+1,idx_start+i-1:idx_start+i+1) = Q(i-1:i+1,i-1:i+1);
end
assignin('base','Q_elastic',Q)


end