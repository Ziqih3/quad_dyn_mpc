classdef Boundary

    properties
        %Boundary
        x_max = 10;
        y_max = 10;
        z_max = 10;
        roll_max = deg2rad(90);
        pitch_max = deg2rad(70);
        yaw_max = deg2rad(90);
        vx_max = 2.0;
        vy_max = 2.0;
        vz_max = 2.0;

        wx_max = pi/2;
        wy_max = pi/2;
        wz_max = pi/2;

        x_min = -10;
        y_min = -10;
        z_min = -10;
        roll_min = -deg2rad(90);
        pitch_min = -deg2rad(70);
        yaw_min = -deg2rad(90);
        vx_min = -2.0;
        vy_min = -2.0;
        vz_min = -2.0;
        wx_min = -pi/2;
        wy_min = -pi/2;
        wz_min = -pi/2;

        % u_min = 1.0975e+05;
        % u_max = 6200000;

        T_min = 0;
        T_max = 200;
        mx_min = -0.3;
        mx_max = 0.3;
        my_min = -0.3;
        my_max = 0.3;
        mz_min = -0.1;
        mz_max = 0.1;
    end

    methods
        function obj = Boundary()
            % Constructor for the limits class
            
        end
    end

end