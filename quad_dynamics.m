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

    tilt_angle = X(10);

    roll_dot_last = X(11);
    pitch_dot_last = X(12);
    yaw_dot_last = X(13);

    thrust_last= X(14);
    
    Thrust = U(1);
    Tilt_speed = U(2);
    rpy_MPC = U(3:5);
    
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
    last_rpy_dot = [roll_dot_last; pitch_dot_last; yaw_dot_last];

    % x(1:3) update equation
    F_aero = CalcAeroForce(V);
    accel = CalcAccel(F_aero,Thrust,rpy,tilt_angle);
    accel = accel + g;

    % x(7:9) update equation
    rpy_sp = rpy_MPC + [0;0;yaw];
    body_rate_sp = diag([2,2,2]) * (rpy_sp - rpy);
    torque = diag([2,2,2]) * (body_rate_sp - rpy_dot) + diag([1,1,1]) * (last_rpy_dot - rpy_dot);

    
    Xdot = SX.sym('Xdot', 14);
    Xdot(1:3) = accel;
    Xdot(4:6) = [roll_dot;pitch_dot;yaw_dot];
    Xdot(7:9) = I_inv*(torque);
    Xdot(10) = Tilt_speed;
    Xdot(11:13) = [roll_dot_last;pitch_dot_last;yaw_dot_last];
    Xdot(14) = thrust_last;

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
%         V = V';
%         a = atand(V(3)/V(1)) + obj.aoi;
aoi = deg2rad(10);
WingSurfaceArea = 0.44;
a = (i/2*log((i+V(3)/V(1)/(i-V(3)/V(1)))) + aoi);

b = 0;
%         R_BW = [cosd(a) * cosd(b),    sind(b),     sind(a)*cosd(b);
%                     -sind(b) * cosd(a),   cosd(b),     -sind(a)*sind(b);
%                     -sind(a)         ,   0  ,        cosd(a)];
R_BW = [cos(a) * cos(b),    sin(b),     sin(a)*cos(b);
    -sin(b) * cos(a),   cos(b),     -sin(a)*sin(b);
    -sin(a)         ,   0  ,        cos(a)];

R_WB = R_BW.';
%       disp((V' * V))
AirDensity = 1.229;          % in kg/m^3
q_bar = (V' * V) * AirDensity / 2;
c_y = 0;
c_z = 0.35 + 0.11 * a;
c_d = 0.01 + 0.2 * a * a;
drag = q_bar * WingSurfaceArea * c_d;
lateral = q_bar * WingSurfaceArea * c_y;
lift = q_bar * WingSurfaceArea * c_z;
%       disp(size(q_bar))
force = R_WB * [-drag; lateral;-lift];
end


