function [Aeq, beq] = getAbeq_bspline1(n, p, start_cond, end_cond, init_type)
% Note: you don't need continuity constraints because the nature of bspline
% garentees you the continuity if you reuse p control points for each segment
% n: n+1 control points
% p: degree of bspline
% note we only need to optimize middle n+1-2*p control points
% the first & last p control points are determined by start & end
% conditions, so start & end conditions dimension are also p x 1
% so in this case, degree = order+1
% first calculate knots
knots = deboor_knot(p,n,2); % we need each curve seg to be clamped at both ends
dknots = knots(2:max(size(knots))-1);
ddknots = dknots(2:max(size(dknots))-1);
dp = p-1;
ddp = dp-1;
dn = n-1;
ddn = dn-1;
C = getC(n,p,knots);
dC = getC(dn,dp,dknots);
ddC = getC(ddn,ddp,ddknots);
assignin('base','C',C)
assignin('base','dC',dC)
assignin('base','ddC',ddC)
%#####################################################
% STEP 2.1 p,v,a constraint in start & end
beq_start_end = [start_cond;end_cond];
Aeq_start = zeros(p, n+1);
Aeq_end = zeros(p, n+1);

C_set = zeros((p-1)*n,n+1);
knots_temp = deboor_knot(p,n,2);
C_current = eye(n+1,n+1);
Bs = zeros(p,p);
Bs(1,1) = 1;
Be = zeros(p,p);
Be(1,p) = 1;
for i = 1:p-1
    p_temp = p-i+1;
    n_temp = n-i+1;
    C_temp = getC(n_temp,p_temp,knots_temp);
    C_current = C_temp * C_current;
    indices_row = 1+(i-1)*n:size(C_current,1)+(i-1)*n;
    indices_coln = 1:size(C_current,2);
    C_set(indices_row,indices_coln) = C_current;
    knots_temp = knots_temp(2:max(size(knots_temp))-1);
    Bs(i+1,1:p) = C_current(1,1:p);
    Be(i+1,1:p) = C_current(n_temp,n-p+2:n+1);
end

Aeq_start(1:p,1:p) = Bs;
Aeq_end(1:p,n-p+2:n+1) = Be;
% Why this is wrong
% because we need all highier orders(v,a,j,...) to change to 0
% continuously,but if we don't set them to 0(set constraints for them),
% they will 
% if init_type == 0
%     Aeq_start = Aeq_start(1,1:end);
%     Aeq_end = Aeq_end(1,1:end);
%     beq_start_end = [start_cond(1);end_cond(1)];
% elseif init_type == 1
%     Aeq_end = Aeq_end(1,1:end);
%     beq_start_end = [start_cond;end_cond(1)];
% end

% Experiment: pos,vel constraint only for start & end
c_num_start = 4; % equal constraint num
c_num_end = 1;
Aeq_start = Aeq_start(1:c_num_start,:);
if c_num_end > 0 
    Aeq_end = Aeq_end(1:c_num_end,:);
    beq_start_end = [start_cond(1:c_num_start);end_cond(1:c_num_end)];
else
    Aeq_end = [];
    beq_start_end = start_cond(1:c_num_start);
end

Aeq = [Aeq_start;Aeq_end];
beq = [beq_start_end];    

end