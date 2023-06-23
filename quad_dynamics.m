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

    tilt_angle = X(10);

    roll_dot_last = X(11);
    pitch_dot_last = X(12);
    yaw_dot_last = X(13);
    
    Thrust = U(1);
    Tilt_speed = U(2);
    rpy_MPC = U(3:5);
    


    %%%% constant for state update
    Mass = 7.4270;
    g = 9.81;
    I_inv = diag([10.685,5.7465,4.6678]);

%%%%%%%%%%%%%%%      state update
    V = [vx;vy;vz];
    rpy = [roll; pitch; yaw];
    rpy_dot = [roll_dot; pitch_dot; yaw_dot];
    last_rpy_dot = [roll_dot_last; pitch_dot_last; yaw_dot_last];

    % x(1:3) update equation
    F_aero = CalcAeroForce(V,rpy);
    accel = CalcAccel(F_aero,Thrust,rpy,tilt_angle)+[0;0;g];


    % x(7:9) update equation
    body_rate_sp = diag([5,5,0])* (rpy_MPC - rpy);
    Torque = diag([5,5,0]) * (body_rate_sp - rpy_dot) + diag([2,2,0]) * (last_rpy_dot - rpy_dot);

    Xdot = [accel;roll_dot;pitch_dot;yaw_dot;I_inv*(Torque);Tilt_speed;roll_dot_last;pitch_dot_last;yaw_dot_last];
%     Xdot(1:3) = accel;
%     Xdot(4:6) = [roll_dot;pitch_dot;yaw_dot];
%     Xdot(7:9) = I_inv*(Torque);
%     Xdot(10) = Tilt_speed;
%     Xdot(11:13) = [roll_dot_last;pitch_dot_last;yaw_dot_last];

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

        RBR = [sin(tilt);
                   0;
              -cos(tilt)];

        T_Body = RBR .* T;

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
V_air = R_BW*V_b;

a_deg = a*180/pi;

q_bar = (V_air'*V_air)*AirDensity / 2;
c_y = 0;
c_z = 0.35 + 0.11 * a_deg;
c_d = 0.03 + 0.2 * a_deg * a_deg;
drag = q_bar*WingSurfaceArea * c_d;
lateral = q_bar*WingSurfaceArea * c_y;
lift = q_bar*WingSurfaceArea *  c_z;

air_speed_norm = norm(V_air);
coeff = aero_blending(air_speed_norm);

drag = coeff * drag;
lateral = coeff * lateral;
lift = coeff * lift;
force = R_WB * [-drag; lateral;-lift];
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




