function xbg = generateSimpleStatebound(N, bound, lower_upper, n_states)

    switch lower_upper
    case 1
        xbg(1:n_states:n_states*(N+1), 1) = bound.x_min;
        xbg(2:n_states:n_states*(N+1), 1) = bound.y_min;
        xbg(3:n_states:n_states*(N+1), 1) = bound.z_min;
        xbg(4:n_states:n_states*(N+1), 1) = bound.roll_min;
        xbg(5:n_states:n_states*(N+1), 1) = bound.pitch_min;
        xbg(6:n_states:n_states*(N+1), 1) = bound.yaw_min;
        xbg(7:n_states:n_states*(N+1), 1) = bound.vx_min;
        xbg(8:n_states:n_states*(N+1), 1) = bound.vy_min;
        xbg(9:n_states:n_states*(N+1), 1) = bound.vz_min;


    case 2
        xbg(1:n_states:n_states*(N+1), 1) = bound.x_max;
        xbg(2:n_states:n_states*(N+1), 1) = bound.y_max;
        xbg(3:n_states:n_states*(N+1), 1) = bound.z_max;
        xbg(4:n_states:n_states*(N+1), 1) = bound.roll_max;
        xbg(5:n_states:n_states*(N+1), 1) = bound.pitch_max;
        xbg(6:n_states:n_states*(N+1), 1) = bound.yaw_max;
        xbg(7:n_states:n_states*(N+1), 1) = bound.vx_max;
        xbg(8:n_states:n_states*(N+1), 1) = bound.vy_max;
        xbg(9:n_states:n_states*(N+1), 1) = bound.vz_max;


    end

end