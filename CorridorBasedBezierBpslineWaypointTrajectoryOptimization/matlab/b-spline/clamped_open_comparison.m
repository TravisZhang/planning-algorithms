% clamped open comparison
% conclusion: ns = n-p+1, there are ns-2(p-1) segs in the middle that are
% the same for open & clamped
close all
clear
c_x = [0 10 20.5 30   40.5 50 60 70 80 90 100]';
c_y = [0 -4 1    6.5  8    10 6  5  10  0  -2]';
pts_knot_list = open('pts_knot.mat');

c_x = [-10 0 10 20.5 30   40.5 50 60 70 80 90 100 110]';
c_y = [-2  0 -4 1    6.5  8    10 6  5  10  0  -2 5]';
pts_knot_list = open('pts_knot1.mat');

c_pts = [c_x c_y];
degree = 4;
dense = 80;
[pts,theta_set,kappa_set,dkappa_set,d_c_pts,dd_c_pts,ddd_c_pts] = get_deboor_attributes(c_pts,degree,dense,1);
[pts1,theta_set1,kappa_set1,dkappa_set1,d_c_pts1,dd_c_pts1,ddd_c_pts1] = get_deboor_attributes(c_pts,degree,dense,2);

figure
plot(pts(:,1),pts(:,2),'.-')
hold on
plot(pts1(:,1),pts1(:,2),'.-')
hold on
plot(c_pts(:,1),c_pts(:,2),'o')
hold on
plot(pts_knot_list.pts_knot_clamped(:,1),pts_knot_list.pts_knot_clamped(:,2),'*')
hold on
plot(pts_knot_list.pts_knot_open(:,1),pts_knot_list.pts_knot_open(:,2),'sq')
axis equal
title('clamped open comparison');