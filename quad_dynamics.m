
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

    
    I_inv = diag([10.685,5.7465,4.6678]);

%%%%%%%%%%%%%%%      state update
    V = [vx;vy;vz];
    rpy = [roll; pitch; yaw];
    rpy_dot = [roll_dot; pitch_dot; yaw_dot];

    % x(1:3) update equation
    F_aero = CalcAeroForce(V);
    
    tilt_angle = 0;
    accel = CalcAccel(F_aero,Thrust,rpy_MPC,tilt_angle);
    accel = accel + g;

    % x(7:9) update equation
    torque = CalcTorque(Thrust,tilt_angle);

    Xdot = SX.sym('Xdot', 9);
    Xdot(1:3) = accel;
    Xdot(4:6) = [roll_dot;pitch_dot;yaw_dot];
    Xdot(7:9) = I_inv*(torque);
    %Xdot(10) = Tilt_speed;
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

Va_i = R_BW*V;
air_speed_norm = norm(Va_i);
q_bar = (Va_i' * Va_i) * AirDensity / 2;

blending_air_speed = 22.194;
Transition_air_speed = 27.7425;

alpha = atan(V(3)/V(1)); %angle attack
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


if air_speed_norm < blending_air_speed
    drag = 0;
    lateral = 0;
    lift = 0;
elseif air_speed_norm >= blending_air_speed && air_speed_norm < Transition_air_speed
    ratio = (air_speed_norm - blending_air_speed) / (Transition_air_speed - blending_air_speed);
    drag = drag_cal * ratio;
    lateral = lateral_cal  * ratio;
    lift = lift_cal * ratio;
elseif air_speed_norm >= Transition_air_speed
    drag = drag_cal;
    lateral = lateral_cal;
    lift = lift_cal;
end

force = [-drag;lateral;-lift];

end

function Torque = CalcTorque(thrust,tilt_angle) %assuming same thrust and tilt for 4 rotor
%ThrustConstant = 1.08105e-4 / 5;    
%TorqueConstant = (1.08105e-4 / 5) * 0.05;

r1 = [0.4;0.4;0];
r2 = [0.4;-0.4;0];
r3 = [-0.4;0.4;0];
r4 = [-0.4;-0.4;0];

T_body = thrust/4*[sin(tilt_angle);0;cos(tilt_angle)];

Torque = cross(r1,T_body)+cross(r2,T_body)+cross(r3,T_body)+cross(r4,T_body);


end

function R_BW = BodyToWind(V)

a = atan(V(3)/V(1)); %angle attack
b = asin(V(2) / norm(V));% angle silp

R_BW = [cos(a) * cos(b),    sin(b),     sin(a)*cos(b);
        -sin(b) * cos(a),   cos(b),     -sin(a)*sin(b);
        -sin(a)         ,   0  ,        cos(a)];
end

function R_RB = RotorToBody(tilt_angle)
R_RB = [cos(tilt_angle), 0, -sin(tilt_angle);
        0,               1,            0;
        sin(tilt_angle), 0, cos(tilt_angle)];

end
