function Xdot = quad_simple(X, U)

    % Use NED frame, 9 states, 4 inputs for standard quad model x, y, z,
    % phi, theta, psi, vx, vy, vz, wx, wy, wz
    addpath('Users/rapstar/Desktop/quad_dyn_mpc/casadi-matlabR2014a-v3.5.5')
    import casadi.*

    x = X(1);
    y = X(2);
    z = X(3);
    
    roll = X(4);
    pitch = X(5);
    yaw = X(6);

    vx = X(7);
    vy = X(8);
    vz = X(9);
    
    Thrust = U(1);
    omega = U(2:4);
    
    R_i2b = GetRotationMatrix(roll, pitch, yaw);
    R_b2i = R_i2b';
    
    nu = [1         sin(roll)*tan(pitch)        cos(roll)*tan(pitch); 
          0         cos(roll)                   -sin(roll);   
          0         sin(roll)/cos(pitch)        cos(roll)/cos(pitch)];
    
    Mass = 7.4270;
    g = 9.81;
%     J = [0.1673   -0.0000         0
%         -0.0000    0.1673         0
%             0         0    0.2879];

    Xdot = SX.sym('Xdot', 9);
    Xdot(1:3) = [vx; vy; vz];
    Xdot(4:6) = nu*omega;
    Xdot(7:9) = Mass*[0.0; 0.0; g] - R_b2i*[0.0; 0.0; Thrust];

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