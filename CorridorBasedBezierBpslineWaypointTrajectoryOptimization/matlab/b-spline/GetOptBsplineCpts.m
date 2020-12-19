function [c_pts_opt,xy_limits] = GetOptBsplineCpts(c_pts,degree,start_cond_xy,end_cond_xy,v_max,a_max,init_type)
% n: n+1 control points
% p: degree of bspline
% Note we need to add 2*p-2 control points to make up for start_cond &
% end_cond determined control points, to fill in middle constraints
n = max(size(c_pts))-1;
p = degree;
corridor_size = 2; % (m)
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
for axis = 1:2
    if init_type == 0 % pos constraint only
        start_cond = zeros(degree,1);
        start_cond(1) = c_pts(1,axis);
        end_cond = zeros(degree,1);
        end_cond(1) = c_pts(n+1,axis); 
    elseif init_type == 1 % full start & end pos
        start_cond = zeros(degree,1);
        start_cond(:) = start_cond_xy(1:degree,axis);
        end_cond = zeros(degree,1);
        end_cond(1) = c_pts(n+1,axis); 
    else % full start & full end
        start_cond = zeros(degree,1);
        start_cond(:) = start_cond_xy(1:degree,axis);
        end_cond = zeros(degree,1);
        end_cond(:) = end_cond_xy(1:degree,axis);
    end
    corridor_range = [-xy_limits(:,axis)+c_pts(:,axis) xy_limits(:,axis)+c_pts(:,axis)];
    axis_name = num2str(axis);
    corridor_name = ['corridor_range_' axis_name];
    assignin('base',corridor_name,corridor_range)
    
    [Aeq, beq] = getAbeq_bspline1(n, p, start_cond, end_cond, init_type);
    assignin('base','Aeq',Aeq)
    assignin('base','beq',beq)
    [Aieq, bieq] = getAbieq_bspline1(n, p, corridor_range, v_max, a_max);
    assignin('base','Aieq',Aieq)
    assignin('base','bieq',bieq)
    [Q] = getQ_bspline(n, p);
    assignin('base','Q',Q)
        
    f = zeros(size(Q,1),1);
    c_pts_opt(:,axis) = quadprog(Q,f,Aieq, bieq, Aeq, beq);
end

end

