% Linear MPC
close all
% init condition
px = 0;
vx = 0;
ax = 0;
jx = 0;
py = 0;
vy = 0;
ay = 0;
jy = 0;
pz = 8;
vz = 0;
az = 0;
jz = 0;
K = 20;
dt = 0.2;
log = [0 px vx ax 0 py vy ay 0 pz vz az 0];
w1 = 1;
w2 = 1;
w3 = 1;
w4 = 1;
w5 = 100; % weight for upper bound slack
w6 = 100; % weight for lower bound slack
vxy_max = 6;
axy_max = 3;
jxy_max = 3;
vz_max = 6;
vz_min = -1;
az_max = 3;
az_min = -1;
jz_max = 2;
jz_min = -2;
total_time = 16;
w_rad = 0.5;
% main process
type = 1; % 0 for no constaints, 1 for hard constraints, 2 for soft constraints
cnt = 1;
for t = 0.2:0.2:total_time
    cnt = cnt+1;
    for i = 1:3
        if i == 1
            p=px;v=vx;a=ax;j=jx;v_max=vxy_max;v_min=-vxy_max;a_max=axy_max;a_min=-axy_max;j_max=jxy_max;j_min=-jxy_max;
        elseif i == 2
            p=py;v=vy;a=ay;j=jy;v_max=vxy_max;v_min=-vxy_max;a_max=axy_max;a_min=-axy_max;j_max=jxy_max;j_min=-jxy_max;
        else
            p=pz;v=vz;a=az;j=jz;v_max=vz_max;v_min=vz_min;a_max=az_max;a_min=az_min;j_max=jz_max;j_min=jz_min;
        end
        % get prediciton matrices
        [Tp,Tv,Ta,Tj,Bp,Bv,Ba,Bj] = getPredictionMaxtrix(K,dt,p,v,a,j);
        % spiral
        p_target = zeros(K,1);
        for k = 1:1:K
            tt = t+k*0.2;
            if(tt > total_time)
                tt = total_time;
            end
            if i == 1
                p_target(k) = 0.5*tt*1.25*cos(w_rad*tt);
            elseif i == 2
                p_target(k) = 0.5*tt*1.25*sin(w_rad*tt);
            else
                p_target(k) = (8-0.5*tt);
            end
        end
        % construct & solve qp problem
        if type == 0
            H = w1*(Tp'*Tp)+w2*(Tv'*Tv)+w3*(Ta'*Ta)+w4*eye(K);
            F = 2*(w1*((Bp-p_target)')*Tp+w2*(Bv')*Tv+w3*(Ba')*Ta);
            J = quadprog(H,F);
        elseif type == 1
            % hard constraints
            H = w1*(Tp'*Tp)+w2*(Tv'*Tv)+w3*(Ta'*Ta)+w4*eye(K);
            F = 2*(w1*((Bp-p_target)')*Tp+w2*(Bv')*Tv+w3*(Ba')*Ta);
            A = [Tv;-Tv;Ta;-Ta;Tj;-Tj];
            b = [ones(K,1)*v_max-Bv;ones(K,1)*(-v_min)+Bv;ones(K,1)*a_max-Ba;ones(K,1)*(-a_min)+Ba;ones(K,1)*j_max-Bj;ones(K,1)*(-j_min)+Bj];
            J = quadprog(H,F,A,b);
        else
            % soft constraints
            H = blkdiag(w1*(Tp'*Tp)+w2*(Tv'*Tv)+w3*(Ta'*Ta)+w4*eye(K),w5*eye(K),w6*eye(K));
            F = [2*(w1*((Bp-p_target)')*Tp+w2*(Bv')*Tv+w3*(Ba')*Ta) zeros(1,K) zeros(1,K)];
            A = [Tv eye(K) zeros(K);-Tv zeros(K) eye(K);Ta eye(K) zeros(K);-Ta zeros(K) eye(K);Tj eye(K) zeros(K);-Tj zeros(K) eye(K)];
            b = [ones(K,1)*v_max-Bv;ones(K,1)*(-v_min)+Bv;ones(K,1)*a_max-Ba;ones(K,1)*(-a_min)+Ba;ones(K,1)*j_max-Bj;ones(K,1)*(-j_min)+Bj];
            J = quadprog(H,F,A,b);
        end

        % apply control to model
        j = J(1);
        p = p+v*dt+1/2*a*dt^2+1/6*j*dt^3;
        v = v+a*dt+1/2*j*dt^2;
        a = a+j*dt;
        % log data
        log(cnt,2+(i-1)*4:1+i*4) = [p v a j];
        if i == 1
            px=p;vx=v;ax=a;jx=j;
        elseif i == 2
            py=p;vy=v;ay=a;jy=j;
        else
            pz=p;vz=v;az=a;jz=j;
        end
    end
    log(cnt,1) = t;
end
% plot data
figure
plot(log(:,1),log(:,2),'.-')
hold on
plot(log(:,1),log(:,3),'.-')
hold on
plot(log(:,1),log(:,4),'.-')
hold on
plot(log(:,1),log(:,5),'.-')
title('x')
legend('p','v','a','j')
figure
plot(log(:,1),log(:,6),'.-')
hold on
plot(log(:,1),log(:,7),'.-')
hold on
plot(log(:,1),log(:,8),'.-')
hold on
plot(log(:,1),log(:,9),'.-')
title('y')
legend('p','v','a','j')
figure
plot(log(:,1),log(:,10),'.-')
hold on
plot(log(:,1),log(:,11),'.-')
hold on
plot(log(:,1),log(:,12),'.-')
hold on
plot(log(:,1),log(:,13),'.-')
title('z')
legend('p','v','a','j')
figure
% plot3(log(1:30,2),log(1:30,6),log(1:30,10))
plot3(log(:,2),log(:,6),log(:,10))
% plot(log(:,2),log(:,2))
