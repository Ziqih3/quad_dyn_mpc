close all

horizon = N;
plot_time = size(xx1, 3);
start_time = 1;


titles = ["vx", "vy", "vz", "Roll", "Pitch", "Yaw", "Roll_dot","pitch_dot", "yaw_dot"...
    , "tilt_angle","Roll_dot_last","pitch_dot_last", "yaw_dot_last"];


figure
for k = 1:12
    subplot(3, 4, k)
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

control_titles = ["Thrust","tilt_speed", "Roll", "Pitch", "Yaw"];
figure
for i = 1:5
    subplot(1, 5, i)
    plot(u_cl(:,i))
    title(control_titles(i))    
end
