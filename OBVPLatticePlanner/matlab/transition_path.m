close all
clear
% transition path calculation
% btwn two paths A & B
% generate path by radius
degree = 4;
dense = 3;
R = 600;
total_dist = 60;
total_num = 30;
offset_x = 50;
offset_y = 10;
total_theta = atan(total_dist/R);
delta_theta = total_theta/total_num;
x_set = zeros(30,1);
y_set = x_set;
for i = 1:30
    theta = i*delta_theta;
    x_set(i) = R*sin(theta);
    y_set(i) = R*(1-cos(theta));
end
% get path A
c_pts0 = [x_set y_set];
[pts0,theta_set0,kappa_set0,dkappa_set0,d_c_pts0,dd_c_pts0,ddd_c_pts0] = get_deboor_attributes(c_pts0,degree,dense,0);
% get path B
c_pts1 = [x_set+offset_x y_set+offset_y];
[pts1,theta_set1,kappa_set1,dkappa_set1,d_c_pts1,dd_c_pts1,ddd_c_pts1] = get_deboor_attributes(c_pts1,degree,dense,0);
total_pt_num = max(size(pts0));

figure
plot(pts0(:,1),pts0(:,2),'-','LineWidth',2)
hold on
plot(pts1(:,1),pts1(:,2),'-','LineWidth',2)
title('start path & target path comparison')
legend('start path A','target path B')
axis equal
%% generate transition path by y = f(x)
% first get start & end points
start_idx = 15;
end_idx = 50;
xs = pts0(start_idx,1);
ys = pts0(start_idx,2);
thetas = theta_set0(start_idx);
ks = kappa_set0(start_idx);
dks = dkappa_set0(start_idx);
xf = pts1(end_idx,1);
yf = pts1(end_idx,2);
thetaf = theta_set1(end_idx);
kf = kappa_set1(end_idx);
dkf = dkappa_set1(end_idx);
% select a region around end point to use for polyfit
idx_range = 10;
range_start = end_idx-idx_range;
range_end = end_idx+idx_range;
if range_end > size(pts1,1)
    range_end = size(pts1,1);
    range_start = range_end - 2*idx_range;
elseif range_start < 1
    range_start = 1;
    range_end = range_start + 2*idx_range; 
end 
polyfit_pts = pts1(range_start:range_end,:);
polyfit_pts_global = pts1(range_start:range_end,:);
% transform end point & polyfit pts to start point coord frame
[xfl,yfl] = global2local(xs,ys,thetas,xf,yf);
thetafl = thetaf - thetas;
for i = 1:size(polyfit_pts,1)
    xg = polyfit_pts(i,1);
    yg = polyfit_pts(i,2);
    [xl,yl] = global2local(xs,ys,thetas,xg,yg);
    polyfit_pts(i,:) = [xl,yl];
end
% calculate transition traj
py0 = 0;
vy0 = 0;
ay0 = ks;
pyf = yfl;
vyf = atan(thetafl);
ayf = kf*(1+vyf^2)^(1.5);
start_end_states = [py0 vy0 ay0 pyf vyf ayf];
[optimal_T, optimal_states] = solve_obvp_x(py0,vy0,ay0,polyfit_pts);
% transform transition pathpoints from local to global
%     t_set = linspace(0,xfl,size(optimal_states,1));
delta_time = 0.05;
t_size = floor(optimal_T/delta_time)+1;
t_set = linspace(0,optimal_T,t_size);
trans_path = zeros(size(optimal_states,1),2);
for i = 1:size(optimal_states,1)
    xl = t_set(i);
    yl = optimal_states(i,1);
    [xg,yg] = local2global(xs,ys,thetas,xl,yl);
    trans_path(i,:) = [xg,yg];
end
% get real end point idx
real_end_idx = 1;
min_dist = sqrt((pts1(1,1)-xg)^2+(pts1(1,2)-yg)^2);
for i = 2:size(pts1,1)
    dist = sqrt((pts1(i,1)-xg)^2+(pts1(i,2)-yg)^2);
    if dist < min_dist
        real_end_idx = i;
        min_dist = dist;
    end
end
% calculate theta kappa dkappa for transition path
theta_tran_set = atan(optimal_states(:,2));
theta_tran_set = theta_tran_set + thetas;
kappa_tran_set = zeros(size(optimal_states,1),1);
for i = 1:size(optimal_states,1)
    dydx = optimal_states(i,2);
    d2ydx2 = optimal_states(i,3);
    kappa_tran_set(i) = compute_curvature_dydx(dydx,d2ydx2);
end

figure
plot(pts0(:,1),pts0(:,2),'-','LineWidth',2)
hold on
plot(pts1(:,1),pts1(:,2),'-','LineWidth',2)
hold on
plot(trans_path(:,1),trans_path(:,2),'-','LineWidth',2)
hold on
plot(polyfit_pts_global(:,1),polyfit_pts_global(:,2),'.-','LineWidth',2)
hold on
plot(trans_path(end,1),trans_path(end,2),'*','LineWidth',2)
title('start path & target path & transition path comparison')
legend('start path A','target path B','transition path','end points','end point')
axis equal
%% functions
function [xl,yl] = global2local(dx,dy,dtheta,xg,yg)
    % B: local, A: global
    ATB = [cos(dtheta) -sin(dtheta) dx;sin(dtheta) cos(dtheta) dy;0 0 1];
    Ap = [xg yg 1].';
    Bp = ATB^-1 * Ap;
    xl = Bp(1);
    yl = Bp(2);
end

function [xg,yg] = local2global(dx,dy,dtheta,xl,yl)
    % B: local, A: global
    ATB = [cos(dtheta) -sin(dtheta) dx;sin(dtheta) cos(dtheta) dy;0 0 1];
    Bp = [xl yl 1].';
    Ap = ATB * Bp;
    xg = Ap(1);
    yg = Ap(2);
end