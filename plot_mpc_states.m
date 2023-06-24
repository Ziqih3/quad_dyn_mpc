close all

horizon = N;
plot_time = size(xx1, 3);
start_time = 1;


titles = ["vx", "vy", "vz", "Roll", "Pitch", "Yaw", "Roll_dot","pitch_dot", "yaw_dot"...
    , "tilt_angle_1","tilt_angle_2","tilt_angle_3","tilt_angle_4","Roll_dot_last","pitch_dot_last"];


figure
for k = 1:15
    subplot(3, 5, k)
    plot(xx(k, :), 'r', 'LineWidth',3)
    hold on
    for i = start_time:plot_time
        t = i:i+horizon;
        z = xx1(:, k, i);
        plot(t, z)
        hold on
    end
    title(titles(k))
end

control_titles = ["Thrust","tilt_speed_1", "tilt_speed_2","tilt_speed_3","tilt_speed_4","roll", "Pitch", "Yaw"];
figure
for i = 1:8
    subplot(2, 4, i)
    plot(u_cl(:,i))
    title(control_titles(i))    
end
