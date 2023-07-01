function Xdot = quad_dynamics(X, U)
    % Use NED frame, 12 states, 4 inputs for standard quad model x, y, z,
    % phi, theta, psi, vx, vy, vz, wx, wy, wz
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

    tilt_angle_1 = X(10);
    tilt_angle_2 = X(11);
    tilt_angle_3 = X(12);
    tilt_angle_4 = X(13);
    tilt_angle = [tilt_angle_1;tilt_angle_2;tilt_angle_3;tilt_angle_4];

    roll_dot_last = X(14);
    pitch_dot_last = X(15);
    yaw_dot_last = X(16);

    Thrust_1 = U(1);
    Thrust_2 = U(2);
    Thrust_3 = U(3);
    Thrust_4 = U(4);
    Thrust = [Thrust_1;Thrust_2;Thrust_3;Thrust_4];

    Tilt_speed_1 = U(5);
    Tilt_speed_2 = U(6);
    Tilt_speed_3 = U(7);
    Tilt_speed_4 = U(8);
    Tilt_speed = [Tilt_speed_1;Tilt_speed_2;Tilt_speed_3;Tilt_speed_4];

    rpy_MPC = U(9:11);
    


    %%%% constant for state update
    Mass = 7.4270;
    g = 9.81;
    I_inv = diag([10.685,5.7465,4.6678]);
    position = [[0.866;0.5;0],[0.866;-0.5;0],[-0.866;-0.5;0],[-0.866;0.5;0]];
%%%%%%%%%%%%%%%      state update
    V = [vx;vy;vz];
    rpy = [roll; pitch; yaw];
    rpy_dot = [roll_dot; pitch_dot; yaw_dot];
    last_rpy_dot = [roll_dot_last; pitch_dot_last; yaw_dot_last];

    % x(1:3) update equation
    F_aero = CalcAeroForce(V,rpy);
    accel = CalcAccel(F_aero,Thrust,rpy,tilt_angle)+[0;0;g];


    % x(7:9) update equation
    %body_rate_sp = diag([5,5,0])* (rpy_MPC - rpy);
    %Torque = diag([5,5,0]) * (body_rate_sp - rpy_dot) + diag([2,2,0]) * (last_rpy_dot - rpy_dot);
    Torque = I_inv*CalcTorque(position,Thrust,tilt_angle);
    Xdot = [accel;roll_dot;pitch_dot;yaw_dot;I_inv*(Torque);Tilt_speed;roll_dot_last;pitch_dot_last;yaw_dot_last];
%     Xdot(1:3) = accel;
%     Xdot(4:6) = [roll_dot;pitch_dot;yaw_dot];
%     Xdot(7:9) = I_inv*(Torque);
%     Xdot(10) = Tilt_speed;
%     Xdot(11:13) = [roll_dot_last;pitch_dot_last;yaw_dot_last];

end


function Torque = CalcTorque(position,Thrust,tilt)
Torque = zeros(3, 1);

RBR_1 = [sin(tilt(1));
    0;
    -cos(tilt(1))];
RBR_2 = [sin(tilt(2));
    0;
    -cos(tilt(2))];
RBR_3 = [sin(tilt(3));
    0;
    -cos(tilt(3))];
RBR_4 = [sin(tilt(4));
    0;
    -cos(tilt(4))];
RBR = [RBR_1,RBR_2,RBR_3,RBR_4];

for i = 1:4
    r = position(:,i);
    T = RBR(:,i)*Thrust(i);
    Torque = Torque + cross(r,T);
end
end

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

function accel = CalcAccel(F_a,T,rpy,tilt)
Mass = 7.4270;

RBR_1 = [sin(tilt(1));
    0;
    -cos(tilt(1))];
RBR_2 = [sin(tilt(2));
    0;
    -cos(tilt(2))];
RBR_3 = [sin(tilt(3));
    0;
    -cos(tilt(3))];
RBR_4 = [sin(tilt(4));
    0;
    -cos(tilt(4))];

T_Body = RBR_1*T(1) + RBR_2*T(2) + RBR_3*T(3) + RBR_4*T(4);

F_tot = F_a + T_Body;

RNI = GetRotationMatrix(rpy(1),rpy(2),rpy(3));
rbi = RNI'; % Body to inertial matrix
F_Inertial = rbi * F_tot;
accel = F_Inertial / Mass;
end

function coeff = aero_blending(air_speed_norm)
% steady fw speed  =  27.7425;
blending_air_speed = 10;
Transition_air_speed = 23;

x = -10:.001:100;

dist = 2;
f1 = @(x) 0;
f2 = @(x) 0.5*(x-blending_air_speed)/(Transition_air_speed-blending_air_speed);
foo = blend(f1, f2, 5, dist);

f3 = @(x) foo(x);
f4 = 1;
foo = blend(f3, f4, 20, dist);

coeff = foo(air_speed_norm);

end
%% piecewise blending
function foo = blend(f1,f2,location,distance)

if nargin < 4
    distance = 0;
end
if nargin < 3 || isempty(location)
    location = 0;
end
validateattributes(location, {'numeric','DimVar'}, {},...
    'blend', 'blending center location', 3);
validateattributes(distance, {'numeric','DimVar'}, {'nonnegative'},...
    'blend', 'blending distance', 4);
if isnumeric(f1)
    f1 = @(~) f1;
end
if isnumeric(f2)
    f2 = @(~) f2;
end
blf = @(x) tanh((x-location)./distance)/2;
foo = @(x) (1/2 - blf(x)).*f1(x) + (1/2 + blf(x)).*f2(x);
end

function force = CalcAeroForce(V,rpy)
R_i2b = GetRotationMatrix(rpy(1),rpy(2),rpy(3));
V_b = R_i2b*V;

a = atan(V_b(3)/V_b(1)); %angle attack;
b = 0;
AirDensity = 1.229;          % in kg/m^3
WingSurfaceArea = 0.44;

R_BW = [cos(a) * cos(b),    sin(b),     sin(a)*cos(b);
    -sin(b) * cos(a),   cos(b),     -sin(a)*sin(b);
    -sin(a)         ,   0  ,        cos(a)];
R_WB = R_BW';

V_air = R_BW * V_b;
a_deg = a*180/pi;

q_bar = (V_air'*V_air) * AirDensity / 2;
c_y = 0;
c_z = 0.35 + 0.11 * a_deg;
c_d = 0.03 + 0.2 * a_deg * a_deg;

drag = q_bar * WingSurfaceArea * c_d;
lateral = q_bar*WingSurfaceArea * c_y;
lift = q_bar*WingSurfaceArea *  c_z;

air_speed_norm = norm(V_air);
coeff = aero_blending(air_speed_norm);

drag = coeff * drag;
lateral = coeff * lateral;
lift = coeff * lift;

force = R_WB * [-drag; lateral;-lift];
force = [0;0;0];
end






