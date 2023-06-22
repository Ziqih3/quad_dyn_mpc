% first casadi test for mpc fpr mobile robots
clear all
close all
clc

import casadi.*

dt = 0.02; %[s]
N = 5; % prediction horizon

% bound = Boundary();

% Velocity
V_x = SX.sym('V_x');
V_y = SX.sym('V_y');
V_z = SX.sym('V_z');
V = [V_x; V_y; V_z];

% Attitude
Roll = SX.sym('Roll');
Pitch = SX.sym('Pitch');
Yaw = SX.sym('Yaw');
rpy = [Roll; Pitch; Yaw];

% Body rate setpoint parametrization
Roll_dot = SX.sym('Roll_dot');
Pitch_dot = SX.sym('Pitch_dot');
Yaw_dot = SX.sym('Yaw_dot');
rpy_dot = [Roll_dot; Pitch_dot; Yaw_dot];

% Tilt angle parametrization
tilt_angle = SX.sym('tilt_angle');

% Last body rate parametrization
Roll_dot_last = SX.sym('Roll_dot_last');
Pitch_dot_last = SX.sym('Pitch_dot_last');
Yaw_dot_last = SX.sym('Yaw_dot_last');
last_rpy_dot = [Roll_dot_last; Pitch_dot_last; Yaw_dot_last];

% 13x1
states = [V; rpy; rpy_dot; tilt_angle; last_rpy_dot];
n_states = length(states);

Roll_MPC = SX.sym('Roll_MPC');
Pitch_MPC = SX.sym('Pitch_MPC');
Yaw_MPC = SX.sym('Yaw_MPC');
rpy_MPC = [Roll_MPC; Pitch_MPC; Yaw_MPC];

% Tilt speed parametrization
tilt_speed = SX.sym('tilt_speed');
% Thrust
Thrust = SX.sym('Thrust');
% Control input
controls= [Thrust; tilt_speed; rpy_MPC];


n_controls = length(controls);

nvar = n_states + n_controls;

rhs = quad_dynamics(states, controls); % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
%P = SX.sym('P',n_states + N*(n_states+n_controls));
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial state and the reference state)

X = SX.sym('X',n_states,(N+1));
% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q_v = diag([20 10 50]);
Q_rpy = diag([10 20 10]);
Q_rpy_dot = diag([10 10 10]);
Q_u = diag([0.0025 1 10 10 10]);
q_ref = [5, 0.1, 0.1]';

st  = X(:,1); % initial state
g = [g;st-P(1:n_states)]; % initial condition constraints
objective_function = 0;

for k = 1:N
    st = X(:,k);  con = U(:,k);
    R_I2B = GetRotationMatrix(st(4),st(5),st(6));
    Velocity_body = R_I2B*st(1:3);

    obj_stateinput = (st(1:3)-P(14:16))' * Q_v * (st(1:3)-P(14:16))...
        +(st(4:6))' * Q_rpy * (st(4:6))...
        +(st(7:9))' * Q_rpy_dot * (st(7:9))...
        + con' * Q_u * con ;
    obj_soft_tilt = exp(-0.332*Velocity_body(1)*(st(10)*180/pi)+1.8*(st(10)*180/pi)+(-0.477)*Velocity_body(1)-2.303);
    objective_function = objective_function +  obj_stateinput + obj_soft_tilt;
   
    st_next = X(:,k+1);
    k1 = f(st, con);   % new
    k2 = f(st + dt/2*k1, con); % new
    k3 = f(st + dt/2*k2, con); % new
    k4 = f(st + dt*k3, con); % new
    st_next_RK4=st +dt/6*(k1 +2*k2 +2*k3 +k4); % new      
    g = [g;st_next-st_next_RK4]; % compute constraints % new
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_controls*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 60000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
% equality constraints

args.lbg(1:n_states*(N+1)) = -1e-20;
args.ubg(1:n_states*(N+1)) = 1e-20;

% input and states constraints
args.lbx(1:n_states:n_states*(N+1),1) = -30; args.ubx(1:n_states:n_states*(N+1),1) = 30; %V in m/s
args.lbx(2:n_states:n_states*(N+1),1) = -30; args.ubx(2:n_states:n_states*(N+1),1) = 30;
args.lbx(3:n_states:n_states*(N+1),1) = -10; args.ubx(3:n_states:n_states*(N+1),1) = 10;
args.lbx(4:n_states:n_states*(N+1),1) = -pi/4; args.ubx(4:n_states:n_states*(N+1),1) = pi/4; %rpy in radian
args.lbx(5:n_states:n_states*(N+1),1) = -pi/4;   args.ubx(5:n_states:n_states*(N+1),1) = pi/6;
args.lbx(6:n_states:n_states*(N+1),1) = -inf; args.ubx(6:n_states:n_states*(N+1),1) = inf;
args.lbx(7:n_states:n_states*(N+1),1) = -pi; args.ubx(7:n_states:n_states*(N+1),1) = pi; %rpy_dot in radian
args.lbx(8:n_states:n_states*(N+1),1) = -pi; args.ubx(8:n_states:n_states*(N+1),1) = pi;
args.lbx(9:n_states:n_states*(N+1),1) = -pi; args.ubx(9:n_states:n_states*(N+1),1) = pi;
args.lbx(10:n_states:n_states*(N+1),1) = -0.122; args.ubx(10:n_states:n_states*(N+1),1) = pi/2; % tilt in degree
args.lbx(11:n_states:n_states*(N+1),1) = -pi; args.ubx(11:n_states:n_states*(N+1),1) = pi;
args.lbx(12:n_states:n_states*(N+1),1) = -pi; args.ubx(12:n_states:n_states*(N+1),1) = pi;
args.lbx(13:n_states:n_states*(N+1),1) = -pi; args.ubx(13:n_states:n_states*(N+1),1) = pi;

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = 0; args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = 90; %thrust in Newton calculated from max pitch angle and weight
args.lbx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = -inf; args.ubx(n_states*(N+1)+2:n_controls:n_states*(N+1)+n_controls*N,1) = inf;
args.lbx(n_states*(N+1)+3:n_controls:n_states*(N+1)+n_controls*N,1) = -pi/4; args.ubx(n_states*(N+1)+3:n_controls:n_states*(N+1)+n_controls*N,1) = pi/4; %rpy_MPC
args.lbx(n_states*(N+1)+4:n_controls:n_states*(N+1)+n_controls*N,1) = -pi/4; args.ubx(n_states*(N+1)+4:n_controls:n_states*(N+1)+n_controls*N,1) = pi/6;
args.lbx(n_states*(N+1)+5:n_controls:n_states*(N+1)+n_controls*N,1) = -pi/4; args.ubx(n_states*(N+1)+5:n_controls:n_states*(N+1)+n_controls*N,1) = pi/4;
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP



% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = zeros(13, 1); x0(1) = 0; x0(2) = 0; x0(3) = 0;% initial state

xs = zeros(13, 1); xs(1) = 0; xs(2) = 0; xs(3) = -5; % goal state

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u_trim = [7.4270*9.81; 0.0; 0.0; 0.0;0.0];

u0 = repmat(u_trim,1,N)'; % control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 15; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
tic
while(norm((x0(1:3)-xs(1:3)),2) > 5e-2 && mpciter < sim_tim / dt)
%while(mpciter < sim_tim / dt)
    %current_time = mpciter * dt;
    args.p   = [x0;xs]; % set the values of the parameters vector
    %----------------------------------------------------------------------
%     args.p(1:n_states) = x0; % initial condition of the robot posture
%     for k = 1:N %new - set the reference to track
%         t_predict = current_time + (k-1)*dt; % predicted time instant
%         z_ref = -0.15*t_predict-2.0; 
%         if z_ref <= -3 % the trajectory end is reached
%             z_ref = -3; 
%         end
%         args.p(nvar*k-3:nvar*k+8) = [0.0, 0.0, z_ref, 0.0, 0.0, 0.0, 0.0, 0.0, -0.15, 0.0, 0.0, 0.0];
%         args.p(nvar*k+9:nvar*k+12) = u_trim;
%     end
    %---------------------------------------------------------------------- 
    % initial value of the optimization variables
    args.x0  = [reshape(X0',n_states*(N+1),1);reshape(u0',n_controls*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    solver.stats().return_status
    u = reshape(full(sol.x(n_states*(N+1)+1:end))',n_controls,N)'; % get controls only from the solution
    xx1(:,1:n_states,mpciter+1)= reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(dt, t0, x0, u,f);
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:n_states*(N+1)))',n_states,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
   
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter
    ss_error = norm((x0(3)-xs(3)),2)
    mpciter = mpciter + 1;
end
toc


%Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam)
plot_mpc_states

%%
%R_i2b
function Rot_BI = GetRotationMatrix(roll, pitch, yaw)

    s_ph = sin(roll);
    s_th = sin(pitch);
    s_ps = sin(yaw);
    c_ph = cos(roll);
    c_th = cos(pitch);
    c_ps = cos(yaw);

    Rot_BI = [ c_th * c_ps                      ,       c_th * s_ps                      ,          -s_th;
               s_ph * s_th * c_ps - c_ph * s_ps ,       s_ph * s_th * s_ps + c_ph * c_ps ,          s_ph * c_th;
               c_ph * s_th * c_ps + s_ph * s_ps ,       c_ph * s_th * s_ps - s_ph * c_ps ,          c_ph * c_th  ];

end
