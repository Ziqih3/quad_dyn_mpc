function xbg = generateStatebound(N, bound, lower_upper)

    switch lower_upper
    case 1
        xbg(1:12:12*(N+1), 1) = bound.x_min;
        xbg(2:12:12*(N+1), 1) = bound.y_min;
        xbg(3:12:12*(N+1), 1) = bound.z_min;
        xbg(4:12:12*(N+1), 1) = bound.roll_min;
        xbg(5:12:12*(N+1), 1) = bound.pitch_min;
        xbg(6:12:12*(N+1), 1) = bound.yaw_min;
        xbg(7:12:12*(N+1), 1) = bound.vx_min;
        xbg(8:12:12*(N+1), 1) = bound.vy_min;
        xbg(9:12:12*(N+1), 1) = bound.vz_min;
        xbg(10:12:12*(N+1), 1) = bound.wx_min;
        xbg(11:12:12*(N+1), 1) = bound.wy_min;
        xbg(12:12:12*(N+1), 1) = bound.wz_min;

    case 2
        xbg(1:12:12*(N+1), 1) = bound.x_max;
        xbg(2:12:12*(N+1), 1) = bound.y_max;
        xbg(3:12:12*(N+1), 1) = bound.z_max;
        xbg(4:12:12*(N+1), 1) = bound.roll_max;
        xbg(5:12:12*(N+1), 1) = bound.pitch_max;
        xbg(6:12:12*(N+1), 1) = bound.yaw_max;
        xbg(7:12:12*(N+1), 1) = bound.vx_max;
        xbg(8:12:12*(N+1), 1) = bound.vy_max;
        xbg(9:12:12*(N+1), 1) = bound.vz_max;
        xbg(10:12:12*(N+1), 1) = bound.wx_max;
        xbg(11:12:12*(N+1), 1) = bound.wy_max;
        xbg(12:12:12*(N+1), 1) = bound.wz_max;


    end

end