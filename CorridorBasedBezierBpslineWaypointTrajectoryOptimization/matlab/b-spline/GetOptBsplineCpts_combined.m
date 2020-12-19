function [c_pts_opt,xy_limits] = GetOptBsplineCpts_combined(c_pts,degree,start_cond_xy,end_cond_xy,v_max,a_max,init_type)
% n: n+1 control points
% p: degree of bspline
% Note we need to add 2*p-2 control points to make up for start_cond &
% end_cond determined control points, to fill in middle constraints
n = max(size(c_pts))-1;
p = degree;
corridor_size = 4; % (m)
xy_limits = zeros(n+1,2);
for i = 1:n
    for j = 1:2
        xy_limits(i,j) = abs(c_pts(i+1,j)-c_pts(i,j));
    end
end
xy_limits(n+1,:) = xy_limits(n,:);
assignin('base','xy_limits',xy_limits)
% we choose max xy_limit element as limits
xy_limit_max = max(xy_limits,[],1);
xy_limits(:,1) = ones(n+1,1)*xy_limit_max(1);
xy_limits(:,2) = ones(n+1,1)*xy_limit_max(2);
xy_limits(:,1) = ones(n+1,1)*corridor_size;
xy_limits(:,2) = ones(n+1,1)*corridor_size;
% last point must be reached(wrong because last/first p points are affected, not
% just one)
% if init_type == 0
% %     xy_limits(1:p,:) = zeros(p,2);
% %     xy_limits(end-p+1:end,:) = zeros(p,2);
%     xy_limits(1,:) = zeros(1,2);
%     xy_limits(end,:) = zeros(1,2);
% elseif init_type == 1
% %     xy_limits(end-p+1:end,:) = zeros(p,2);
%     xy_limits(end,:) = zeros(1,2);
% end
assignin('base','xy_limits_real',xy_limits)

c_pts_opt = zeros(n+1,2);

if init_type == 0 % pos constraint only
    start_cond = zeros(degree,2);
    start_cond(1,:) = c_pts(1,:);
    end_cond = zeros(degree,2);
    end_cond(1,:) = c_pts(n+1,:); 
elseif init_type == 1 % full start & end pos
    start_cond = zeros(degree,2);
    start_cond(1:degree,:) = start_cond_xy(1:degree,:);
    end_cond = zeros(degree,2);
    end_cond(1,:) = c_pts(n+1,:); 
else % full start & full end
    start_cond = zeros(degree,2);
    start_cond(1:degree,:) = start_cond_xy(1:degree,:);
    end_cond = zeros(degree,2);
    end_cond(1:degree,:) = end_cond_xy(1:degree,:);
end
% min x, max x, min y, max y
corridor_range = [-xy_limits(:,1)+c_pts(:,1) xy_limits(:,1)+c_pts(:,1) -xy_limits(:,2)+c_pts(:,2) xy_limits(:,2)+c_pts(:,2)];
corridor_name = ['corridor_range'];
assignin('base',corridor_name,corridor_range)

[Aeq, beq] = getAbeq_bspline_combined(n, p, start_cond, end_cond, init_type);
assignin('base','Aeq',Aeq)
assignin('base','beq',beq)
[Aieq, bieq] = getAbieq_bspline_combined(n, p, corridor_range, v_max, a_max);
assignin('base','Aieq',Aieq)
assignin('base','bieq',bieq)
[Q] = getQ_bspline_combined(n, p);
assignin('base','Q',Q)

f = zeros(size(Q,1),1);
c_pts_opt_temp = quadprog(Q,f,Aieq, bieq, Aeq, beq);
c_pts_opt(1:n+1,1) = c_pts_opt_temp(1:n+1);
c_pts_opt(1:n+1,2) = c_pts_opt_temp(n+2:2*(n+1));

end

