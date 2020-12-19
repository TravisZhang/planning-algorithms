close all
% clear
files = dir('*.csv');
% plot 3 consecutive paths from sim
c_pts0 = [currentglobal0.x currentglobal0.y];
c_pts1 = [currentglobal1.x currentglobal1.y];
c_pts2 = [currentglobal2.x currentglobal2.y];

dense = 3;
degree = 4;
[pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(c_pts0,degree,dense,1);
[pts1,theta_set1,kappa_set1,dkappa_set1,d_c_pts1,dd_c_pts1,ddd_c_pts1] = get_deboor_attributes(c_pts1,degree,dense,1);
[pts2,theta_set2,kappa_set2,dkappa_set2,d_c_pts2,dd_c_pts2,ddd_c_pts2] = get_deboor_attributes(c_pts2,degree,dense,1);

% 2,3 ... paths for start idx in path 1
start_idx_2 = find_pt_idx(pts,pts1(1,:));
start_idx_3 = find_pt_idx(pts,pts2(1,:));
start_idx_2 = start_idx_2-1;
start_idx_3 = start_idx_3-1;

% plot 3 paths in one
figure
plot(pts(:,1),pts(:,2),'.-')
hold on
plot(pts1(:,1),pts1(:,2),'.-')
hold on
plot(pts2(:,1),pts2(:,2),'.-')
hold on
plot(c_pts0(:,1),c_pts0(:,2),'o')
hold on
plot(c_pts1(:,1),c_pts1(:,2),'+')
hold on
plot(c_pts2(:,1),c_pts2(:,2),'sq','MarkerSize',10)
axis equal
title('comparison of 3 paths, start clamp for 1, no clamp for 2/3');
legend('path0','path1','path2','cpts0','cpts1','cpts2')
% plot theta/kappa/dkappa for 3 paths
idx_set = linspace(1,size(pts,1),size(pts,1));
figure
plot(theta_set*180/pi, '.-')
hold on
plot(idx_set+start_idx_2,theta_set1*180/pi, '+-')
hold on
plot(idx_set+start_idx_3,theta_set2*180/pi, 'o-')
title('theta set in degree');
figure
plot(kappa_set, '.-')
hold on
plot(idx_set+start_idx_2,kappa_set1, '+-')
hold on
plot(idx_set+start_idx_3,kappa_set2, 'o-')
title('kappa set')
figure
plot(dkappa_set, '.-')
hold on
plot(idx_set+start_idx_2,dkappa_set1, '+-')
hold on
plot(idx_set+start_idx_3,dkappa_set2, 'o-')
title('dkappa set')
%% test opt comparison
c_pts_origin0 = [cptsorigin0.x cptsorigin0.y];
c_pts_origin1 = [cptsorigin1.x cptsorigin1.y];
c_pts_origin2 = [cptsorigin2.x cptsorigin2.y];
c_pts_opt_0 = [cptsopt0.x cptsopt0.y];
c_pts_opt_1 = [cptsopt1.x cptsopt1.y];
c_pts_opt_2 = [cptsopt2.x cptsopt2.y];
[output] = OptMain(c_pts_origin0,degree,dense);
[output1] = OptMain(c_pts_origin1,degree,dense);
[output2] = OptMain(c_pts_origin2,degree,dense);

figure
plot(output.c_pts_opt(:,1),output.c_pts_opt(:,2),'.-')
hold on
plot(c_pts_opt_0(:,1),c_pts_opt_0(:,2),'.-')
title('comparison for cptsopt set 0')
legend('matlab','real')
figure
plot(output1.c_pts_opt(:,1),output1.c_pts_opt(:,2),'.-')
hold on
plot(c_pts_opt_1(:,1),c_pts_opt_1(:,2),'.-')
title('comparison for cptsopt set 1')
legend('matlab','real')
figure
plot(output2.c_pts_opt(:,1),output2.c_pts_opt(:,2),'.-')
hold on
plot(c_pts_opt_2(:,1),c_pts_opt_2(:,2),'.-')
title('comparison for cptsopt set 2')
legend('matlab','real')

% compare real opt with matlab opt
[pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(c_pts_opt_0,degree,dense,1);
[pts1,theta_set1,kappa_set1,dkappa_set1,d_c_pts1,dd_c_pts1,ddd_c_pts1] = get_deboor_attributes(output.c_pts_opt,degree,dense,1);

figure
subplot(3,1,1)
plot(theta_set*180/pi, '.-')
hold on
plot(theta_set1*180/pi, '+-')
title('theta set in degree');
legend('real','matlab')
subplot(3,1,2)
plot(kappa_set, '.-')
hold on
plot(kappa_set1, '+-')
title('kappa set')
legend('real','matlab')
subplot(3,1,3)
plot(dkappa_set, '.-')
hold on
plot(dkappa_set1, '+-')
title('dkappa set')
legend('real','matlab')

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