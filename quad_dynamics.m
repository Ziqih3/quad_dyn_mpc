
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
    

    R_i2b = GetRotationMatrix(roll, pitch, yaw);
    R_b2i = R_i2b';
    
    nu = [1         sin(roll)*tan(pitch)        cos(roll)*tan(pitch); 
          0         cos(roll)                   -sin(roll);   
          0         sin(roll)/cos(pitch)        cos(roll)/cos(pitch)];

    %%%% constant for state update
    Mass = 7.4270;
    g = 9.81;
    J = [0.1673   -0.0000         0
        -0.0000    0.1673         0
            0         0    0.2879];
    I_inv = diag([10.685,5.7465,4.6678]);

%%%%%%%%%%%%%%%      state update
    V = [vx;vy;vz];
    rpy = [roll; pitch; yaw];
    rpy_dot = [roll_dot; pitch_dot; yaw_dot];

    % x(1:3) update equation
    F_aero = CalcAeroForce(V);
    tilt_angle = 0;

    accel = CalcAccel(F_aero,Thrust,rpy,tilt_angle);
    accel = accel + g;

    % x(7:9) update equation
    rpy_sp = rpy_MPC + [0;0;yaw];
    body_rate_sp = diag([2,2,2]) * (rpy_sp - rpy);
    torque = diag([2,2,2]) * (body_rate_sp - rpy_dot) + diag([1,1,1]) * (last_rpy_dot - rpy_dot);


    Xdot = SX.sym('Xdot', 13);
    Xdot(1:3) = accel;
    Xdot(4:6) = [roll_dot;pitch_dot;yaw_dot];
    Xdot(7:9) = I_inv*(torque);
    %Xdot(10) = Tilt_speed;
    Xdot(10:12) = [roll_dot;pitch_dot;yaw_dot];
    Xdot(13) = Thrust;

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

function accel = CalcAccel(F_a,T,rpy,tilt)
        
        Mass = 7.4270;

        RBR = [sin(tilt);
                   0;
              -cos(tilt)];
  
        T_Body = RBR .* T;

        F_tot = F_a + T_Body;
        angles = [rpy(3) rpy(2) rpy(1)];
        
        cang = cos(angles);
        sang = sin(angles);

        RNI = [cang(2).*cang(1)                            , cang(2).*sang(1)                            , -sang(2);
               sang(3).*sang(2).*cang(1) - cang(3).*sang(1), sang(3).*sang(2).*sang(1) + cang(3).*cang(1), sang(3).*cang(2);
               cang(3).*sang(2).*cang(1) + sang(3).*sang(1), cang(3).*sang(2).*sang(1) - sang(3).*cang(1), cang(3).*cang(2)];

        rbi = RNI'; % Body to inertial matrix
        F_Inertial = rbi * F_tot;

        accel = F_Inertial / Mass;
end


function force = CalcAeroForce(V)

WingSurfaceArea = 0.44;
AirDensity = 1.229;         
R_BW = BodyToWind(V);

V_Air = R_BW*V;

q_bar = (V_Air' * V_Air) * AirDensity / 2;

c_y = 0;
c_z = 0.35 + 0.11 * a;
c_d = 0.01 + 0.2 * a * a;
drag = q_bar * WingSurfaceArea * c_d;
lateral = q_bar * WingSurfaceArea * c_y;
lift = q_bar * WingSurfaceArea * c_z;

R_WB = R_BW';              % Wind to body
force = R_WB * [-drag; lateral;-lift];
end

function Torque = CalcAerodynamicTorque(V) %TODO
a = atan(V(3)/V(1)); %angle attack

R_BW = BodyToWind(V);
V_Air = R_BW*V;
q_bar = (V_Air' * V_Air) * physics.AirDensity / 2;

%dimensionless angular rate
Wingspan = 2; % in meters
MeanChord = 0.22; % in meters
WingSurfaceArea = 0.44;

%Torque = [L,M,N] body frame
cl = get_cl(a);
cm = get_cm(a);

roll_torque = q_bar * WingSurfaceArea * Wingspan * cl;
pitch_torque = q_bar * WingSurfaceArea * MeanChord * cm;
yaw_torque = q_bar * WingSurfaceArea * Wingspan * 0;

Torque = [roll_torque; pitch_torque; yaw_torque];
end

function R_BW = BodyToWind(V)

a = atan(V(3)/V(1)); %angle attack
b = asin(V(2) / norm(V));% angle silp

R_BW = [cos(a) * cos(b),    sin(b),     sin(a)*cos(b);
        -sin(b) * cos(a),   cos(b),     -sin(a)*sin(b);
        -sin(a)         ,   0  ,        cos(a)];
end

function cl = get_cl(alpha)
    data = load('C_l');
    cl = fixpt_interp1([-8.5:.25:13.75, 14.5,14.75, 15],data.C_l,alpha,sfix(8),2^-3,sfix(16), 2^-14,'Nearest');
end


function cm = get_cm(alpha)
data = load('C_m');
cm = fixpt_interp1([-8.5:.25:13.75, 14.5,14.75, 15],data.C_m,alpha,sfix(8),2^-3,sfix(16), 2^-14,'Floor');
end
