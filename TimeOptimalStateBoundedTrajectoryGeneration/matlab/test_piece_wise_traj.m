%% test piece wise traj
close all
clear

% Ref: On-Line Planning of Time-Optimal, Jerk-Limited Trajectories

%% 2nd order traj
% with 3 phases(acc,cruise,dec to halt)
x0 = -3;
xg = 1.5; % 50 for no solution with adjusted time 8s
        % 2 for below 2.4(zero cruise time)
v0 = -1; % 4 -4 (exceed v limit)
v_max = 3;
a_max = 2;
[dt1,dt2,dt3,states] = get_2nd_order_traj_profile(x0,xg,v0,v_max,a_max);

% we need to unify total time for different dimensions
% choose the longest total time within all dimensions
% and adjust dt1 dt2 dt3 according to new(longer) total time T
% decrease dt1,dt3 both by dt, increase dt2'=T-dt1'-dt3'
T = 8;
[dt1n,dt2n,dt3n,statesn] = get_adjusted_2nd_order_traj_profile(dt1,dt2,dt3,x0,xg,v0,v_max,a_max,T);

% my own derivation of delta time for adjusment(proved to be correct)
enable_derivation = 0;
if enable_derivation
    syms v v0 dt1 dt2 dt3 a_max d T delta
    a_acc = d*a_max;
    a_dec = -d*a_max;
    dt1n = dt1 - delta;
    dt3n = dt3 - delta;
    dt2n = T - dt1n - dt3n;
    v = v0 + a_acc*dt1;
    [dx1,~] = second_order_traj(dt1, 0, v0, a_acc);
    [dx2,~] = second_order_traj(dt2, 0, v, 0);
    [dx3,~] = second_order_traj(dt3, 0, v, a_dec);
    total_dist = dx1+dx2+dx3;
    vn = v0+a_acc*dt1n;
    LHS = 2*vn^2-v0^2+2*d*a_max*vn*dt2n-2*d*a_max*total_dist;
    [coeffs_lhs,delta_set] = coeffs(LHS,delta);
    roots_matlab = roots(coeffs_lhs);
    roots_delta = get_roots(flip(coeffs_lhs));
    roots_delta = simplify(roots_delta);

    % method in paper(appear to be incorrect ??)
    t3 = dt1+dt2+dt3;
    A = T-(t3-dt2);
    delta_example = -A/2+sqrt(A^2/4+(T-t3)*v0/a_max);
end
 

%% 3rd order traj
a0 = -1;
j_max = 3;

% zero cruise profile
[dt_set,states_zero] = zero_cruise_profile(x0,xg,v0,a0,v_max,a_max,j_max);

% time optimal traj
[states3,dt_set3,input3] = get_3rd_order_traj_profile(x0,xg,v0,a0,v_max,a_max,j_max);

% given time traj
[states3n,dt_set3n,input3n] = get_adjusted_3rd_order_traj_profile(x0,xg,v0,a0,v_max,a_max,j_max,T);

% time optimal traj(set speed)
vg = 2;
[states3_setspeed,input3_setspeed] = get_setspeed_3rd_order_traj_profile(x0,v0,vg,a0,v_max,a_max,j_max);

% given time traj(set speed)
[states3_setspeedn,input3_setspeedn] = get_setspeed_adjusted_3rd_order_traj_profile(x0,v0,vg,a0,v_max,a_max,j_max,T);

enable_derivation_3rd_order = 0;
if enable_derivation_3rd_order
    syms a a0 v_max a_max j_max d d0 d1 xg x0
    j_pos = d0*j_max;
    j_neg = -d0*j_max;
    a = d0*a_max;
    dt1 = (a-a0)/j_pos;
    [dv1,~] = second_order_traj(dt1, 0, a0, j_pos);
    [dx1,v1,a1] = third_order_traj(dt1, 0, v0, a0, j_pos);
    dt3 = -a/j_neg;
    [dv3,~] = second_order_traj(dt3, 0, a, j_neg);    
    dt2 = (d*v_max-(v0+dv1+dv3))/a;
    [dx2,v2,a2] = third_order_traj(dt2, 0, v1, a, 0);
    [dx3,v3,a3] = third_order_traj(dt3, 0, v2, a, j_neg);

    j_pos = d1*j_max;
    j_neg = -d1*j_max;
    aa = d1*a_max;
    dt4 = (aa-0)/j_pos;
    [dv4,~] = second_order_traj(dt4, 0, 0, j_pos);
    [dx4,v4,a4] = third_order_traj(dt4, 0, d*v_max, a0, j_pos);
    dt6 = -aa/j_neg;
    [dv6,~] = second_order_traj(dt6, 0, aa, j_neg);    
    dt5 = (0-(d*v_max+dv4+dv6))/aa;
    [dx5,v5,a5] = third_order_traj(dt5, 0, v4, aa, 0);
    [dx6,v6,a6] = third_order_traj(dt6, 0, v5, aa, j_neg);
    
    ds = xg-x0-(dx1+dx2+dx3+dx4+dx5);
    dsdv_max = diff(ds,v_max);
    [coeffs_v_max,v_max_set] = coeffs(dsdv_max,v_max);
    roots_v_max = get_roots(flip(coeffs_v_max));
end

%% plot
% plot_states(states)
% plot_states(statesn)
plot_states_3rd_order(states_zero)
plot_states_3rd_order(states3)
plot_states_3rd_order(states3n)
plot_states_3rd_order(states3_setspeed)
plot_states_3rd_order(states3_setspeedn)

%% functions
function plot_states(states)

figure
plot(states(:,1),states(:,2),'-')
hold on
plot(states(:,1),states(:,3),'-')
hold on
plot(states(:,1),states(:,4),'-')
title(['traj profile, total time: ',num2str(states(end,1)), ' s'])
legend('pos(m)','vel(m/s)','acc(m/s^2)')

end

function plot_states_3rd_order(states)

figure
plot(states(:,1),states(:,2),'-')
hold on
plot(states(:,1),states(:,3),'-')
hold on
plot(states(:,1),states(:,4),'-')
hold on
plot(states(:,1),states(:,5),'-')
title(['traj profile, total time: ',num2str(states(end,1)), ' s'])
legend('pos(m)','vel(m/s)','acc(m/s^2)','jerk(m/s^3)')

end

function [dt1n,dt2n,dt3n,states] = get_adjusted_2nd_order_traj_profile(dt1,dt2,dt3,x0,xg,v0,v_max,a_max,T)

% first judge sign of v_max and a_max
% by computing full stop position xs and compare with xg
[d] = get_acc_vel_sign(x0,xg,v0,a_max);

% compute delta time change for phase 1 & 3
inside_sqr = T^2*a_max^2*d^2 + 2*T*a_max^2*d^2*dt1 - 2*T*a_max^2*d^2*dt3 + 4*T*a_max*d*v0 - a_max^2*d^2*dt1^2 - 6*a_max^2*d^2*dt1*dt3 - 4*dt2*a_max^2*d^2*dt1 + 3*a_max^2*d^2*dt3^2 - 8*a_max*d*dt3*v0 - 4*dt2*a_max*d*v0 + 2*v0^2;
% note: if inside_sqr < 0, then there's no solution !!!!
delta1 = ((inside_sqr)^(1/2) - T*a_max*d + a_max*d*dt1 + a_max*d*dt3)/(2*a_max*d);
delta2 = -((inside_sqr)^(1/2) + T*a_max*d - a_max*d*dt1 - a_max*d*dt3)/(2*a_max*d);
if delta1 <= 0 
    delta = delta2;
elseif delta2 <= 0
    delta = delta1;
else
    delta = min(delta1,delta2);
end
delta_paper = dt1/2 - T/2 + dt3/2 + ((dt1 - T + dt3)^2/4 - (v0*(dt1 - T + dt2 + dt3))/a_max)^(1/2);
dt1n = dt1 - delta;
dt3n = dt3 - delta;
dt2n = T-dt1n-dt3n;
% get new delta time for 3 phases
a_acc = d*a_max;
a_dec = -d*a_max;
% for delta > dt1, then the first phase has to be a dec phase
% dec -> cruise -> dec
if delta > dt1
    ts = (0-v0)/as;
    vn = (xg-xs)/(T-ts); % cruise dist / cruise time
    a_dec = -d*a_max;
    a_acc = a_dec;
    dt1n = (vn-v0)/a_acc;
    dt3n = -vn/a_dec;
    dt2n = T-dt1n-dt3n;
end

% get states
ts = 0.01;
[states] = get_traj_states_2nd_order(dt1n,dt2n,dt3n,x0,v0,a_acc,a_dec,ts);

end


function [dt1,dt2,dt3,states,input,d] = get_2nd_order_traj_profile(x0,xg,v0,v_max,a_max)

% first judge sign of v_max and a_max
% by computing full stop position xs and compare with xg
[d] = get_acc_vel_sign(x0,xg,v0,a_max);

% get delta time for 3 phases
a_acc = d*a_max;
a_dec = -d*a_max;
v = d*v_max;
dt1 = (v-v0)/a_acc;
[dx1,~] = second_order_traj(dt1, 0, v0, a_acc); % relative pos
dt3 = -v/a_dec;
[dx3,~] = second_order_traj(dt3, 0, v, a_dec); % relative pos
dt2 = (xg-(x0+dx1+dx3))/v;
if dt2 <= 0 % there's no cruise phase
    v_abs = sqrt(d*a_max*(xg-x0)+0.5*v0^2);
    v = d*v_abs;
    dt2 = 0;
    dt1 = (v-v0)/a_acc;
    dt3 = -v/a_dec;
end

% get states
ts = 0.01;
[states,input] = get_traj_states_2nd_order(dt1,dt2,dt3,x0,v0,a_acc,a_dec,ts);

end


function [states,input] = get_setspeed_adjusted_3rd_order_traj_profile(x0,v0,vg,a0,v_max,a_max,j_max,T)

% generate vel profile, state=[v,a,j], input=j
[dt1,dt2,dt3,states,input,d] = get_2nd_order_traj_profile(v0,vg,a0,a_max,j_max);
T_total = dt1+dt2+dt3;
if T > T_total
    [dt1,dt2,dt3,states] = get_adjusted_2nd_order_traj_profile(dt1,dt2,dt3,v0,vg,a0,a_max,j_max,T);
end

dt_set = [dt1,dt2,dt3];

% get states(pos profile)
ts = 0.01;
[states,input] = get_traj_states_3rd_order(dt_set,x0,v0,a0,input,ts);

end


function [states,input] = get_setspeed_3rd_order_traj_profile(x0,v0,vg,a0,v_max,a_max,j_max)

% generate vel profile, state=[v,a,j], input=j
[dt1,dt2,dt3,states,input,d] = get_2nd_order_traj_profile(v0,vg,a0,a_max,j_max);

dt_set = [dt1,dt2,dt3];

% get states(pos profile)
ts = 0.01;
[states,input] = get_traj_states_3rd_order(dt_set,x0,v0,a0,input,ts);

end


function [states,dt_set,input] = get_adjusted_3rd_order_traj_profile(x0,xg,v0,a0,v_max,a_max,j_max,T)

T_total = inf;
v3 = v_max;

% if traj too fast, reduce v3
search_method = 1;
if search_method == 0 % line search
    max_iter = 500;
    dt_threshold = 0.01;
    v3_reduce_quantity = 0.01;
    iter = 0;
    result_set = [];
    while abs(T-T_total) > dt_threshold && iter < max_iter
        [states,dt_set,input,v3] = get_3rd_order_traj_profile(x0,xg,v0,a0,v3,a_max,j_max);
        T_total = sum(dt_set);
        if T < T_total
            break
        end
        v3 = v3 - v3_reduce_quantity;
        iter = iter+1;
        result_set = [result_set;[iter,T-T_total,v3]];
    end
else % equal interval search (much faster)
    [states,dt_set,input,v3] = get_3rd_order_traj_profile(x0,xg,v0,a0,v3,a_max,j_max);
    T_total = sum(dt_set);
    result_set = [];
    if T > T_total % else if time optimal traj already too slow at start, then faster traj doesn't exist
        max_iter = 500;
        dt_threshold = 0.01;
        iter = 0;
        epsilon = 0.01; % it's crucial to add interval
        va = 0;
        vb = v_max;
        while abs(T-T_total) > dt_threshold && iter < max_iter
            v3 = min((va+vb)/2,v_max);
            [states,dt_set,input,v3] = get_3rd_order_traj_profile(x0,xg,v0,a0,v3,a_max,j_max);
            T_total = sum(dt_set);
            if T > T_total % too fast, reduce v3
                vb = v3+epsilon;
            else
                va = v3-epsilon;
            end
            iter = iter+1;
            result_set = [result_set;[iter,T-T_total,v3]];
        end
    end
end

assignin('base','result_set_time_opt',result_set)

end


function [states,dt_set,input,v3_final] = get_3rd_order_traj_profile(x0,xg,v0,a0,v_max,a_max,j_max)

ds = inf;
v3 = v_max;

% select correct profile type
% if pos overshoot, reduce v3, else add cruise duration
search_method = 1;
if search_method == 0 % line search
    max_iter = 500;
    ds_threshold = 0.001;
    v3_reduce_quantity = 0.01;
    iter = 0;
    result_set = [];
    while abs(ds) > ds_threshold && iter < max_iter
        [d, ds, dt_set, states, input] = get_acc_vel_sign_3rd_order(x0,xg,v0,a0,v3,a_max,j_max);
        if d > 0 % add cruising duration
            [dt_set,input] = add_cruising_duration(dt_set, input, ds, d*v3);
            % get states(pos profile)
            ts = 0.01;
            [states,input] = get_traj_states_3rd_order(dt_set,x0,v0,a0,input,ts);
            break;
        else % reduce v3
            v3 = v3 - v3_reduce_quantity;
    %         d0 = d_set(1);d1 = d_set(2);d = d_set(3);
    %         v_max_root = ((d*((a0*(a0 - a_max*d0))/(d0*j_max) - (a0 - a_max*d0)^2/(2*d0*j_max) + 1))/(a_max*d0) - (2*a_max*d)/j_max + (d*((a0 - a_max*d0)^2/(2*d0*j_max) + (a_max^2*d0)/(2*j_max) - (a0*(a0 - a_max*d0))/(d0*j_max) - 1))/(a_max*d0) + (d*((d1*a_max^2)/(2*j_max) + (a0*a_max)/j_max))/(a_max*d1))/(d^2/(a_max*d0) - d^2/(a_max*d1));
        end
        iter = iter+1;
        result_set = [result_set;[iter,ds,v3]];
    end
else % equal interval search (much faster)
    max_iter = 500;
    ds_threshold = 0.001;
    iter = 0;
    result_set = [];
    epsilon = 0.01; % it's crucial to add interval
    [d, ds, dt_set, states, input] = get_acc_vel_sign_3rd_order(x0,xg,v0,a0,v3,a_max,j_max);
    if d > 0 % add cruising duration
        [dt_set,input] = add_cruising_duration(dt_set, input, ds, d*v3);
        % get states(pos profile)
        ts = 0.01;
        [states,input] = get_traj_states_3rd_order(dt_set,x0,v0,a0,input,ts);
    else % xs(stop pos) overshoot xg
        va = 0;
        vb = v_max;
        while abs(ds) > ds_threshold && iter < max_iter
            v3 = min((va+vb)/2,v_max);
            [d, ds, dt_set, states, input] = get_acc_vel_sign_3rd_order(x0,xg,v0,a0,v3,a_max,j_max);
            if d < 0 % reduce v3
                vb = v3+epsilon;
            else
                va = v3-epsilon;
            end  
            iter = iter+1;
            result_set = [result_set;[iter,ds,v3]];
        end
    end
end

v3_final = v3;

assignin('base','result_set',result_set)

fprintf(1,'3rd order time optimal, iters: %d, ds: %.3f\n',iter,ds)

end


function [dt_set,input] = add_cruising_duration(dt_set, input, ds, v)
    dt4 = ds / v;
    input = [input(1:3),0,input(4:6)];
    dt_set = [dt_set(1:3),dt4,dt_set(4:6)];
end


function [d] = get_acc_vel_sign(x0,xg,v0,a_max)

ds = xg-x0;
if v0 ~= 0
    dv = 0 - v0;
    as = a_max*sign(dv);
    xs = final_pos_second_order_traj(x0,v0,0,as);
    ds = xg-xs;
end

d = sign(ds);

end

% the paper suggest to just use 2nd order one
function [d, ds, dt_set, states, input] = get_acc_vel_sign_3rd_order(x0,xg,v0,a0,v_max,a_max,j_max)

ds = xg-x0;
dz = 1;
if v0 ~= 0 || a0 ~= 0
    % compute zero-cruise profile
    [dt_set,states,input,d_set] = zero_cruise_profile(x0,xg,v0,a0,v_max,a_max,j_max);
    xs = states(end,2);
    ds = xg-xs;
    dz = d_set(3);
end

d = sign(ds)*dz;

end


function [dt_set,states,input,d_set] = zero_cruise_profile(x0,xg,v0,a0,v_max,a_max,j_max)
% Note:
% 1. We just use 2nd order profile but change the state from [p,v,a] to [v,a,j]
% 2. There are two vel profiles to compute, 
%      a. s0=[v0,a0,aj] to s1 = [v_max, 0, 0]
%      b. s1 to s2 = [0,0,0]
%    The two profiles combine to get zero-cruise profile
% 3. The position is accumulated as mentioned

% [d] = get_acc_vel_sign(x0,xg,v0,a_max); % inaccurate when a0 != 0
% instead we use 3rd order decrease to zero profile
[states,input] = get_setspeed_3rd_order_traj_profile(x0,v0,0,a0,v_max,a_max,j_max);
d = sign(xg-states(end,2));

% 1st vel profile(s0 -> s1)
[dt1a,dt2a,dt3a,statesa,inputa,d0] = get_2nd_order_traj_profile(v0,d*v_max,a0,a_max,j_max);

% 2nd vel profile(s0 -> s1)
[dt1b,dt2b,dt3b,statesb,inputb,d1] = get_2nd_order_traj_profile(d*v_max,0,0,a_max,j_max);

% combine two vel profiles
dt_set = [dt1a,dt2a,dt3a,dt1b,dt2b,dt3b];
input = [inputa,inputb];

% plot vel traj
% statesb(:,1) = statesb(:,1) + statesa(end,1);
% states = [statesa;statesb];
% plot_states(states)

% get states(pos profile)
ts = 0.01;
[states,input] = get_traj_states_3rd_order(dt_set,x0,v0,a0,input,ts);

d_set = [d0,d1,d];
end


function [states,input] = get_traj_states_2nd_order(dt1,dt2,dt3,x0,v0,a_acc,a_dec,ts)

t3 = dt1+dt2+dt3;
t2 = dt1+dt2;
t1 = dt1;
v = v0+dt1*a_acc;
[dx1,~] = second_order_traj(dt1, 0, v0, a_acc);
[dx2,~] = second_order_traj(dt2, 0, v, 0);
x1 = x0+dx1;
x2 = x1+dx2;
states = [];
for t = 0:ts:t3
    if t < t1
        [xt,vt] = second_order_traj(t, x0, v0, a_acc);
        state = [t,xt,vt,a_acc];
    elseif t < t2
        [xt,vt] = second_order_traj(t-t1, x1, v, 0);
        state = [t,xt,vt,0];
    else
        [xt,vt] = second_order_traj(t-t2, x2, v, a_dec);
        state = [t,xt,vt,a_dec];
    end
    states = [states;state];
end

input = [a_acc,0,a_dec];

end


function [states,input] = get_traj_states_3rd_order(dt_set,x0,v0,a0,input,ts)

states = [];
t0 = 0;
last_state = [x0,v0,a0];
for i = 1:max(size(dt_set))
    if dt_set(i) == 0
        continue
    end
    for t = 0:ts:dt_set(i)
        [xt,vt,at] = third_order_traj(t, last_state(1), last_state(2), last_state(3), input(i));
        state = [t+t0,xt,vt,at,input(i)];
        states = [states;state];
    end
    [xt,vt,at] = third_order_traj(dt_set(i), last_state(1), last_state(2), last_state(3), input(i));
    last_state = [xt,vt,at];
    t0 = t0+dt_set(i);
end

end


function [x] = final_pos_second_order_traj(x0,v0,vt,a)

x = (vt^2-v0^2)/2/a+x0;

end


function [xt,vt] = second_order_traj(dt, x0, v0, a)

xt = x0 + v0*dt + 1/2*a*dt^2;
vt = v0 + a*dt;

end

function [xt,vt,at] = third_order_traj(dt, x0, v0, a0, j)

xt = x0 + v0*dt + 1/2*a0*dt^2 + 1/6*j*dt^3;
vt = v0 + a0*dt + 1/2*j*dt^2;
at = a0 + j*dt;

end


function root_set = get_roots(coeffs)

coeffs = coeffs./coeffs(max(size(coeffs))); % we need to normalize coeff(n_order)
n = max(size(coeffs))-1;
CM = sym(diag(ones(n-1,1),-1));
for i = 1:n
    CM(i,n) = -coeffs(i);
end
fprintf(1,'rank of C: %d\n',rank(CM));
root_set = eig(CM); % possible values of xf/T
assignin('base','roots',root_set)

end