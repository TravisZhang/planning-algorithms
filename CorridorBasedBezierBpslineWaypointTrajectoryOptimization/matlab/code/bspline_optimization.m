close all
clear
%% first segment of path to optimize
c_x = [0 10 20.5 30   40.5 50 60 70 80 90 100]';
c_y = [0 -4 1    6.5  8    10 6  5  10  0  -2]';
c_pts = [c_x c_y];
degree = 4;
dense = 80;
[output] = OptMain(c_pts,degree,dense);
%% second segment of path to optimize
% we need to reuse at least 4 points(=degree) to ensure 3rd order
% continuity
c_pts1 = output.c_pts_opt(end-degree+1:end,:); % minimal reuse
c_pts1 = output.c_pts_opt(end-degree-3:end,:); % more reuse
c_x1 = [110 120 130 140 150 160].';
c_y1 = [4   8   2   -2  -6  1].';
c_pts1 = [c_pts1;[c_x1 c_y1]];
% c_pts1 = [c_x1 c_y1];
[output1] = OptMain(c_pts1,degree,dense);
% start_cond_xy = zeros(4,2);
% for i = 1:2
%     start_cond_xy(:,i) = [output.c_pts_opt(end,i);output.d_c_pts_opt(end,i);output.dd_c_pts_opt(end,i);output.ddd_c_pts_opt(1,i)];
% end
% [output1] = OptMain(c_pts1,degree,dense,start_cond_xy);
%% third segment of path to optimize
% we need to reuse at least 4 points(=degree) to ensure 3rd order
% continuity
c_pts2 = output1.c_pts_opt(end-degree+1:end,:); % minimal reuse
c_x2 = [170 180 190 200 210 220].';
c_y2 = [10   2   -4   5  -5  1].';
c_pts2 = [c_pts2;[c_x2 c_y2]];
% c_pts1 = [c_x1 c_y1];
[output2] = OptMain(c_pts2,degree,dense);
% start_cond_xy = zeros(4,2);
% for i = 1:2
%     start_cond_xy(:,i) = [output.c_pts_opt(end,i);output.d_c_pts_opt(end,i);output.dd_c_pts_opt(end,i);output.ddd_c_pts_opt(1,i)];
% end
% [output1] = OptMain(c_pts1,degree,dense,start_cond_xy);
%% comparison of two paths
% figure
% plot(output.pts_opt(:,1),output.pts_opt(:,2),'.-')
% hold on
% plot(output1.pts_opt(:,1),output1.pts_opt(:,2),'.-')
% hold on
% plot(output2.pts_opt(:,1),output2.pts_opt(:,2),'.-')
% axis equal
% title('comparison of two paths');
%% change knot type for path 2/3
[pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(output.c_pts_opt,degree,dense,0);
% [pts1,theta_set1,kappa_set1,dkappa_set1,d_c_pts1,dd_c_pts1,ddd_c_pts1] = get_deboor_attributes(output1.c_pts_opt,degree,dense,1);
c_pts_opt_1 = [output.c_pts_opt(1:end-degree-3-1,:);output1.c_pts_opt];
c_pts_opt_1 = output1.c_pts_opt;
[pts1,theta_set1,kappa_set1,dkappa_set1,d_c_pts1,dd_c_pts1,ddd_c_pts1] = get_deboor_attributes(c_pts_opt_1,degree,dense,1);
[pts2,theta_set2,kappa_set2,dkappa_set2,d_c_pts2,dd_c_pts2,ddd_c_pts2] = get_deboor_attributes(output2.c_pts_opt,degree,dense,1);
% 2,3 ... paths for start idx in path 1
start_idx_2 = find_pt_idx(pts,pts1(1,:));
start_idx_3 = find_pt_idx(pts1,pts2(1,:));
start_idx_2 = start_idx_2-1;
start_idx_3 = start_idx_3+start_idx_2-1;

% plot 3 paths in one
figure
plot(pts(:,1),pts(:,2),'.-')
hold on
plot(pts1(:,1),pts1(:,2),'.-')
hold on
plot(pts2(:,1),pts2(:,2),'.-')
hold on
plot(output.c_pts_opt(:,1),output.c_pts_opt(:,2),'o')
hold on
plot(c_pts_opt_1(:,1),c_pts_opt_1(:,2),'+')
hold on
plot(output2.c_pts_opt(:,1),output2.c_pts_opt(:,2),'sq','MarkerSize',10)
axis equal
title('comparison of 3 paths, start clamp for 1, no clamp for 2/3');
legend('path0','path1','path2','cpts0','cpts1','cpts2')
% plot theta/kappa/dkappa for 3 paths
idx_set = linspace(1,size(pts,1),size(pts,1));
idx_set1 = linspace(1,size(pts1,1),size(pts1,1));
idx_set2 = linspace(1,size(pts2,1),size(pts2,1));

figure
plot(theta_set*180/pi, '.-')
hold on
plot(idx_set1+start_idx_2,theta_set1*180/pi, '+-')
hold on
plot(idx_set2+start_idx_3,theta_set2*180/pi, 'o-')
title('theta set in degree');
figure
plot(kappa_set, '.-')
hold on
plot(idx_set1+start_idx_2,kappa_set1, '+-')
hold on
plot(idx_set2+start_idx_3,kappa_set2, 'o-')
title('kappa set')
figure
plot(dkappa_set, '.-')
hold on
plot(idx_set1+start_idx_2,dkappa_set1, '+-')
hold on
plot(idx_set2+start_idx_3,dkappa_set2, 'o-')
title('dkappa set')

%% using Aeq & start_cond/end_cond to reverse calculate cpts
% this indicates that Aeq is the same for same degree!!!
% conclusion: the whole curve will be different, including c_pts, knots,
% etc. but we only know the first p points(p=degree) according to
% start_cond, so for now it's not possible to generate new path i.e. the same
% as old one, the only way is maybe to rearrange the middle n+1-2*p of c_pts & knots.
Cs = Aeq(1:degree,1:degree);
Ce = Aeq(degree+1:2*degree, end-degree+1:end);
% we choose a point on 1st bspline to be start point
start_idx = 30;
[pos_xy,vel_xy,acc_xy,jerk_xy] = match_point(c_pts,degree,dense,0,start_idx);
c_pts_new = zeros(degree,2);
for i = 1:2
    c_pts_new(1:degree,i) = Cs^(-1)*[pos_xy(i) vel_xy(i) acc_xy(i) jerk_xy(i)].';
end
c_pts_new = [c_pts_new;c_pts(degree+1:end,:)];
[pts3,theta_set3,kappa_set3,dkappa_set3,d_c_pts3,dd_c_pts3,ddd_c_pts3] = get_deboor_attributes(c_pts_new,degree,dense,0);
[pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(c_pts,degree,dense,0);
% optimiza new set
[output3] = OptMain(c_pts_new,degree,dense);
[pts4,theta_set4,kappa_set4,dkappa_set4,d_c_pts4,dd_c_pts4,ddd_c_pts4] = get_deboor_attributes(output3.c_pts_opt,degree,dense,0);
% plot new points compare to previous points
figure
plot(pts(:,1),pts(:,2),'.-')
hold on
plot(c_pts(:,1),c_pts(:,2),'o')
hold on
plot(pts3(:,1),pts3(:,2),'.-')
hold on
plot(c_pts_new(:,1),c_pts_new(:,2),'+')
hold on
plot(pts4(:,1),pts4(:,2),'.-')
hold on
plot(output3.c_pts_opt(:,1),output3.c_pts_opt(:,2),'sq')
hold on
axis equal
title('comparison of new cpts with old cpts')
legend('origin pts','origin c_pts','new pts','new c_pts','opt pts','opt c_pts')


%% functions
function idx = find_pt_idx(path,pt)
idx = 1;
dist = xy_dist(path(1,:),pt);
for i = 1:size(path,1)
    dist_temp = xy_dist(path(i,:),pt);
    if dist_temp < dist
        dist = dist_temp;
        idx = i;
    end
end

end

function dist = xy_dist(pt0,pt1)
dist = sqrt((pt0(1)-pt1(1))^2+(pt0(2)-pt1(2))^2);
end