%***************************************
%Author: Chaoyu Zhang
%Date: 2020-04
%***************************************
%% ���̳�ʼ��
clear all; close all;
x_I=1; y_I=1;           % ���ó�ʼ��
x_G=700; y_G=700;       % ����Ŀ���
Thr=50;                 %����Ŀ�����ֵ
Delta= 30;              % ������չ����
FLAG_REWIRE = 1;        % RRT* ��֦����
rewire_radius = 2*Delta;% RRT*
FLAG_INFORMED = 1;      % Informed RRT* (generate random new point in a Ellipse that its const_dist = current path dist)
flag_path_found = 0;    % informed RRT* can only perform when a path is found
goal_index = 0;         % informed RRT* has to know last point to goal
flag_informed_done = 0; % in case we are not generating random point inside ellipse
max_iter_in = 2000;     % max iter number for inform random
prob_toward_goal = 0.8; % if p > prob_toward_goal select goal as new point
% Explanation:
% the const dist of ellipse = PA+PB, P is a point on curve, A and B are its
% two two focal points
%% ������ʼ��
T.v(1).x = x_I;         % T������Ҫ��������v�ǽڵ㣬�����Ȱ���ʼ����뵽T������
T.v(1).y = y_I; 
T.v(1).xPrev = x_I;     % ��ʼ�ڵ�ĸ��ڵ���Ȼ���䱾��
T.v(1).yPrev = y_I;
T.v(1).dist=0;          %�Ӹ��ڵ㵽�ýڵ�ľ��룬�����ȡŷ�Ͼ���
T.v(1).indPrev = 0;     %
%% ��ʼ������������ҵ����
figure(1);
ImpRgb=imread('newmap.png');
Imp=rgb2gray(ImpRgb);
imshow(Imp)
xL=size(Imp,1);%��ͼx�᳤��
yL=size(Imp,2);%��ͼy�᳤��
hold on
plot(x_I, y_I, 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
plot(x_G, y_G, 'go', 'MarkerSize',10, 'MarkerFaceColor','g');% ��������Ŀ���
count=1;
record_path_length = 0;
for iter = 1:3000
    x_rand=[];
    %Step 1: �ڵ�ͼ���������һ����x_rand
    %��ʾ���ã�x_rand(1),x_rand(2)����ʾ�����в����������
    x_rand = zeros(1,2);
    x_rand(1) = floor(xL * rand()) + 1;
    x_rand(2) = floor(yL * rand()) + 1;
    % Informed RRT* random sample in ellipse
    p_rand = rand();
    if p_rand > prob_toward_goal
        x_rand = [x_G y_G];
    elseif FLAG_INFORMED == 1 && flag_path_found == 1
        iter_in = 0;
        while 1
            x_rand = zeros(1,2);
            x_rand(1) = floor(xL * rand()) + 1;
            x_rand(2) = floor(yL * rand()) + 1;
            last_point = [T.v(goal_index).x T.v(goal_index).y];
            dist_temp = cal_dist(last_point,[x_G y_G]);
            const_dist = T.v(goal_index).dist + dist_temp;
            record_path_length = const_dist;
            a = const_dist/2; % half long axis
            c = cal_dist([x_I y_I],[x_G y_G])/2; % focal of ellipse
            b = sqrt(a^2-c^2); % half short axis
            theta = atan2(y_G-y_I,x_G-x_I);
            coord_local = [cos(theta) -sin(theta) x_I;sin(theta) cos(theta) y_I;0 0 1]^-1 * [x_rand.';1];
            coord_local(1) = coord_local(1) - c;
            if coord_local(1)^2/a^2 + coord_local(2)^2/b^2 <= 1
                break
            end
            if iter_in > max_iter_in
                flag_informed_done = 1;
%                 break;
            end
            iter_in = iter_in+1;
        end
    end
    x_near=[];
    %Step 2: ���������������ҵ�����ڽ���x_near 
    %��ʾ��x_near�Ѿ�����T��
    x_near = [T.v(1).x T.v(1).y];
    min_dist = cal_dist(x_near,x_rand);
    x_near_idx = 1;
    for i = 1:size(T.v,2)
        x_near_temp = [T.v(i).x T.v(i).y];
        dist_temp = cal_dist(x_near_temp,x_rand);
        if dist_temp < min_dist
            x_near = x_near_temp;
            x_near_idx = i;
            min_dist = dist_temp;
        end
    end
    x_new=[];
    %Step 3: ��չ�õ�x_new�ڵ�
    %��ʾ��ע��ʹ����չ����Delta
    theta = atan2(x_rand(2)-x_near(2),x_rand(1)-x_near(1));
    x_new = [Delta*cos(theta) Delta*sin(theta)] + x_near;
    %���ڵ��Ƿ���collision-free
    if ~collisionChecking(x_near,x_new,Imp) 
       continue;
    end
    count=count+1;
    
    % Rewire for RRT*
    % cal path length from start to x_new
    c_path_dist = Delta + T.v(x_near_idx).dist;
    if FLAG_REWIRE == 1 || FLAG_INFORMED == 1
        % find neighbor nodes in rewire_radius
        for i = 1:size(T.v,2)
            if i == x_near_idx
                continue % skip x_near
            end
            x_near_temp = [T.v(i).x T.v(i).y];
            dist_temp = cal_dist(x_near_temp,x_new);
            if dist_temp < rewire_radius
                % rewire if path from x_near_i is shorter than its current
                % path from x_near
                n_path_dist = dist_temp + T.v(i).dist;
                if c_path_dist > n_path_dist
                    x_near_idx = i;
                    c_path_dist = n_path_dist;
                    x_near = [T.v(i).x T.v(i).y];
                end
            end
        end 
    end
    
    %Step 4: ��x_new������T 
    %��ʾ���½ڵ�x_new�ĸ��ڵ���x_near
    current_size = size(T.v,2);
    new_idx = current_size+1;
    T.v(new_idx).x = x_new(1);         % T������Ҫ��������v�ǽڵ�
    T.v(new_idx).y = x_new(2); 
    T.v(new_idx).xPrev = x_near(1);     
    T.v(new_idx).yPrev = x_near(2);
    T.v(new_idx).dist = c_path_dist; %�Ӹ��ڵ㵽�ýڵ�ľ��룬�����ȡŷ�Ͼ���
    T.v(new_idx).indPrev = x_near_idx;     %
    if FLAG_INFORMED == 1 && flag_path_found == 1
        index_list = [goal_index];
        pathIndex = goal_index;
        while 1
            pathIndex = T.v(pathIndex).indPrev;
            index_list = [index_list;pathIndex];
            if pathIndex == 1
                break
            end
        end
        for i = length(index_list)-1:-1:1
            idx = index_list(i);
            idx_p = index_list(i+1);
            dist_temp = cal_dist([T.v(idx).x T.v(idx).y],[T.v(idx_p).x T.v(idx_p).y]);
            T.v(idx).dist = T.v(idx_p).dist + dist_temp;
        end
    end
    %Step 5:����Ƿ񵽴�Ŀ��㸽�� 
    %��ʾ��ע��ʹ��Ŀ�����ֵThr������ǰ�ڵ���յ��ŷʽ����С��Thr����������ǰforѭ��
    fprintf(1,'path_length: %.3f\n',record_path_length);
    goal_dist = cal_dist(x_new, [x_G y_G]);
    if goal_dist < Thr 
        if FLAG_INFORMED == 0
            break
        else
            flag_path_found = 1;
            goal_index = new_idx;
        end
    end
   %Step 6:��x_near��x_new֮���·��������
   %��ʾ 1��ʹ��plot���ƣ���ΪҪ�����ͬһ��ͼ�ϻ����߶Σ�����ÿ��ʹ��plot����Ҫ����hold on����
   %��ʾ 2�����ж��յ���������forѭ��ǰ���ǵð�x_near��x_new֮���·��������
   line_seg = [x_near;x_new];
   plot(line_seg(:,1),line_seg(:,2),'r.-')
   hold on
   pause(0.005); %��ͣ0.1s��ʹ��RRT��չ�������׹۲�
end
%% ·���Ѿ��ҵ��������ѯ
if iter < 2000 || flag_path_found == 1
    path.pos(1).x = x_G; path.pos(1).y = y_G;
    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
    pathIndex = T.v(end).indPrev; % �յ����·��
    j=0;
    while 1
        path.pos(j+3).x = T.v(pathIndex).x;
        path.pos(j+3).y = T.v(pathIndex).y;
        pathIndex = T.v(pathIndex).indPrev;
        if pathIndex == 1
            break
        end
        j=j+1;
    end  % ���յ���ݵ����
    path.pos(end+1).x = x_I; path.pos(end).y = y_I; % ������·��
    for j = 2:length(path.pos)
        plot([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y], 'b', 'Linewidth', 3);
    end
else
    disp('Error, no path found!');
end


