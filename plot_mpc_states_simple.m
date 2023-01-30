close all

horizon = N;
plot_time = size(xx1, 3);
start_time = 1;


titles = ["X", "Y", "Z", "Roll", "Pitch", "Yaw", "Vx", "Vy", "Vz"];


figure
for k = 1:9
    subplot(3, 3, k)
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

control_titles = ["Thrust", "wx", "wy", "wz"];
figure
for i = 1:4
    subplot(2, 2, i)
    plot(u_cl(:,i))
    title(control_titles(i))    
end
