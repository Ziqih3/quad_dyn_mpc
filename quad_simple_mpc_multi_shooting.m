
% first casadi test for mpc fpr mobile robots
clear all
close all
clc

addpath('Users/rapstar/Desktop/quad_dyn_mpc/casadi-matlabR2014a-v3.5.5')
import casadi.*

dt = 0.02; %[s]
N = 50; % prediction horizon

bound = Boundary();

x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');
roll = SX.sym('roll'); pitch = SX.sym('pitch'); yaw = SX.sym('yaw');
vx = SX.sym('vx'); vy = SX.sym('vy'); vz = SX.sym('vz');

states = [x; y; z; roll; pitch; yaw; vx; vy; vz];
n_states = length(states);

thrust = SX.sym('thrust'); wx = SX.sym('wx'); wy = SX.sym('wy'); wz = SX.sym('wz');
controls = [thrust; wx; wy; wz]; 
n_controls = length(controls);

nvar = n_states + n_controls;

rhs = quad_simple(states, controls); % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
%P = SX.sym('P',n_states + N*(n_states+n_controls));
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial state and the reference state)

X = SX.sym('X',n_states,(N+1));
% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = eye(n_states,n_states); Q(1:3,1:3) = diag([1;1;2]); Q(4:6,4:6) = 0.5*eye(3); Q(7:9,7:9) = diag([1;1;2]);% weighing matrices (states)
R = eye(n_controls,n_controls); R(1,1) = 0.5; R(2:4,2:4) = 0.05*eye(3); % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:n_states)]; % initial condition constraints
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(n_states+1:2*n_states))'*Q*(st-P(n_states+1:2*n_states)) + con'*R*con; % calculate obj
%     obj = obj + (st - P(nvar*k-3:nvar*k+8))'*Q*(st - P(nvar*k-3:nvar*k+8)) + ...
%           (con - P(nvar*k+9:nvar*k+12))'*R*(con - P(nvar*k+9:nvar*k+12));
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

args.lbg(1:n_states*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_states*(N+1)) = 0;  % 1e-20   % Equality constraints

args.lbx = generateSimpleStatebound(N, bound, 1, n_states);
args.ubx = generateSimpleStatebound(N, bound, 2, n_states);

args.lbx((n_states*(N+1) + 1):4:(n_states*(N+1) + n_controls*N),1) = bound.T_min;
args.ubx((n_states*(N+1) + 1):4:(n_states*(N+1) + n_controls*N),1) = bound.T_max;
args.lbx((n_states*(N+1) + 2):4:(n_states*(N+1) + n_controls*N),1) = bound.wx_min;
args.ubx((n_states*(N+1) + 2):4:(n_states*(N+1) + n_controls*N),1) = bound.wx_max;
args.lbx((n_states*(N+1) + 3):4:(n_states*(N+1) + n_controls*N),1) = bound.wy_min;
args.ubx((n_states*(N+1) + 3):4:(n_states*(N+1) + n_controls*N),1) = bound.wy_max;
args.lbx((n_states*(N+1) + 4):4:(n_states*(N+1) + n_controls*N),1) = bound.wz_min;
args.ubx((n_states*(N+1) + 4):4:(n_states*(N+1) + n_controls*N),1) = bound.wz_max;
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = zeros(n_states, 1); x0(1) = 0; x0(2) = 0; x0(3) = -2.0;% initial state
xs = zeros(n_states, 1); xs(1) = 1; xs(2) = -2; xs(3) = -4.0; xs(6) = deg2rad(-10);% goal state

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u_trim = [7.4270*9.81; 0.0; 0.0; 0.0];

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
%while(norm((x0(3)-xs(3)),2) + norm((x0(6)-xs(6)),2) > 5e-2 && mpciter < sim_tim / dt)
while(norm((x0-xs),2) > 5e-2 && mpciter < sim_tim / dt)
%while(~all((x0-xs) < 1e-5) && mpciter < sim_tim / dt)
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
    ss_error = norm((x0(3)-xs(3)),2);
    mpciter = mpciter + 1;
end
toc


%Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam)
plot_mpc_states_simple

