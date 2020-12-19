close all
clear

state_num = 6;
graph = ones(state_num,state_num).*(-1);
% plug in % one-step costs for edges with direction
% if cost >= 0, then the action is valid
% else if cost < 0, it means the cost is inf
% Ss = s5, Sg = s6
% state can have multiple robot actions (u), and one action can have
% multiple nature actions (theta)
graph(1,2) = 2;
graph(1,6) = 2;
graph(2,1) = 2;
graph(2,4) = 1;
graph(3,6) = 1;
graph(4,3) = 3;
graph(5,2) = 1;
% set theta edges
graph_theta = zeros(state_num,state_num);
graph_theta(1,2) = 1;
graph_theta(1,6) = 1;

goal_state = 6;
start_state = 5;

% find optimal cost-to-goal G* for every state
% also known as nondeterministic dijkstra
G_set = ones(state_num,1).*(-1); % -1 means inf
G_set(goal_state) = 0;
state_status = zeros(state_num,1); % 0 for unknown, 1 for open set, -1 for closed set
state_status(goal_state) = 1; % put xf in open set
while 1
    x1 = find_state_with_min_g_value(G_set,state_status);
    state_status(x1) = -1;
    if x1 == start_state
        fprintf(1,'reached xs !!!\n')
        break
    end
    x0_set = find_predecessor(graph, x1);
    for i = 1:max(size(x0_set))
        x0 = x0_set(i);
        if state_status(x0) == -1
            continue
        end
        max_G = calculate_max_g_value(graph, graph_theta, x0, G_set);
        if G_set(x0) > max_G || G_set(x0) == -1 && max_G ~= -1
            G_set(x0) = max_G;
            state_status(x0) = 1;
        end
    end
end

function max_G = calculate_max_g_value(graph, graph_theta, x0, G_set)
% find max cost-to-goal for x0
    x1_set = find_successor_theta(graph, graph_theta, x0);
    if isempty(x1_set)
        x1 = graph(x0,
        max_G = 
    end
    max_G = -1;
    for i = 1:max(size(x1_set))
        x1 = x1_set(i);
        G1 = G_set(x1);
        l = graph(x0,x1); % one-step cost
        if l == -1
            
        end
        if G1 == -1
            return % inf is definitely max value
        end
        G_temp = l+G1;
        if G_temp > max_G
            max_G = G_temp;
        end
    end
end

function x0_set = find_predecessor(graph, x1)
    x0_set = [];
    for s = 1:max(size(graph))
        if graph(s,x1) > -1
            x0_set = [x0_set;s];
        end
    end
end

function x1_set = find_successor_theta(graph, graph_theta, x0)
    x1_set = [];
    for s = 1:max(size(graph))
        if graph(x0,s) > -1 && graph_theta(x0,s) == 1
            x1_set = [x1_set;s];
        end
    end
end

function [x] = find_state_with_min_g_value(G_set,state_status)
    min_G = -1;
    x = 1;
    for s = 1:max(size(state_status))
        if state_status(s) ~= 1 %  not in open set
            continue
        end
        G_temp = G_set(s);
        if min_G > G_temp || (min_G == -1 && G_temp > -1)
            min_G = G_temp;
            x = s;
        end
    end
end