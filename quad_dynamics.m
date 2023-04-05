
function Xdot = quad_dynamics(X, U)

    % Use NED frame, 12 states, 4 inputs for standard quad model x, y, z,
    % phi, theta, psi, vx, vy, vz, wx, wy, wz
    addpath('/Users/junyigeng/Downloads/casadi-osx-matlabR2015a-v3.5.5')
    import casadi.*

    vx = X(1);
    vy = X(2);
    vz = X(3);

    roll = X(4);
    pitch = X(5);
    yaw = X(6);

    roll_dot = X(7);
    pitch_dot = X(8);
    yaw_dot = X(9);

    %tilt_angle = X(10);
    
    Thrust = U(1);
    %Tilt_speed = U(2);
    rpy_MPC = U(2:4);
    %%%% constant for state update
    Mass = 7.4270;
    g = 9.81;
    I_inv = diag([10.685,5.7465,4.6678]);

%%%%%%%%%%%%%%%      state update
    V = [vx;vy;vz];
    rpy = [roll; pitch; yaw];
    rpy_dot = [roll_dot; pitch_dot; yaw_dot];

    % x(1:3) update equation
    F_aero = CalcAeroForce(V);
    tilt_angle = 0;
    accel = CalcAccel(F_aero,Thrust,rpy,tilt_angle) + [0;0;g];
    torque = CalcTorque(rpy_MPC,rpy,rpy_dot);
    % x(7:9) update equation

    Xdot = SX.sym('Xdot', 9);
    Xdot(1:3) = accel;
    Xdot(4:6) = rpy_dot;
    Xdot(7:9) = I_inv*(torque);
    %Xdot(10) = Tilt_speed;
end





function accel = CalcAccel(F_a,T,rpy,tilt)
        Mass = 7.4270;
        RBR = [sin(tilt);
                   0;
              -cos(tilt)];
        T_Body = RBR*T;
        F_tot = F_a + T_Body;
        R_i2b = GetRotationMatrix(rpy(1),rpy(2),rpy(3));
        R_b2i = R_i2b'; 
        F_Inertial = R_b2i * F_tot;
        accel = F_Inertial / Mass;
end


function force = CalcAeroForce(V)

WingSurfaceArea = 0.44;
AirDensity = 1.229; 

R_BW = BodyToWind(V);
Va_i = R_BW*V;

air_speed_norm = norm(Va_i);
q_bar = (Va_i' * Va_i) * AirDensity / 2;

blending_air_speed = 22.194;
Transition_air_speed = 27.7425;  

alpha = atan(V(3)/V(1))*180/pi; %angle attack
C_Z0 = 0.35;
C_Za = 0.11;
C_D0 = 0.01;
C_Da = 0.2;

c_z = C_Z0 + C_Za * alpha;
c_x = C_D0 + C_Da * alpha * alpha;
c_y = 0;

drag_cal = q_bar *WingSurfaceArea * c_x;
lateral_cal = q_bar * WingSurfaceArea * c_y;
lift_cal = q_bar * WingSurfaceArea * c_z;

force = R_BW'*[-drag_cal;lateral_cal;-lift_cal];
%force = [0;0;0];

end

function Torque = CalcTorque(rpy_MPC,rpy,rpy_dot) %assuming same thrust and tilt for 4 rotor
%ThrustConstant = 1.08105e-4 / 5;    
%TorqueConstant = (1.08105e-4 / 5) * 0.05;
% r1 = [0.4;0.4;0];
% r2 = [0.4;-0.4;0];
% r3 = [-0.4;0.4;0];
% r4 = [-0.4;-0.4;0];
% T_body = thrust/4*[sin(tilt_angle);0;cos(tilt_angle)];
% Torque = cross(r1,T_body)+cross(r2,T_body)+cross(r3,T_body)+cross(r4,T_body);
rpy_sp = rpy_MPC ;
body_rate_sp = diag([2,2,0])* (rpy_sp - rpy);
Torque = diag([2,2,0]) * (body_rate_sp - rpy_dot);

end

function R_BW = BodyToWind(V)

a = atan(V(3)/V(1)); %angle attack
b = 0;% angle silp

R_BW = [cos(a) * cos(b),    sin(b),     sin(a)*cos(b);
        -sin(b) * cos(a),   cos(b),     -sin(a)*sin(b);
        -sin(a)         ,   0  ,        cos(a)];
end


function R_RB = RotorToBody(tilt_angle)
R_RB = [cos(tilt_angle), 0, -sin(tilt_angle);
        0,               1,            0;
        sin(tilt_angle), 0, cos(tilt_angle)];

end

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