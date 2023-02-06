close all

horizon = N;
plot_time = size(xx1, 3);
start_time = 1;


titles = ["vx", "vy", "vz", "Roll", "Pitch", "Yaw", "Roll_dot","pitch_dot", "yaw_dot"...
    , "Roll_dot_last","pitch_dot_last", "yaw_dot_last",'last_thrust'];
%"tilt_angle",

xx = [xx(1:3,:);rad2deg(xx(4:12,:));xx(13,:)];
figure
for k = 1:12
    subplot(3, 4, k)
    plot(xx(k, :), 'r', 'LineWidth',3)
    hold on
%     for i = start_time:plot_time
%         t = i:i+horizon;
%         z = xx1(:, k, i);
%         plot(t, z)
%         hold on
%     end
    title(titles(k))
end
%subplot()
u_cl = [u_cl(:,1),rad2deg(u_cl(:,2:end))];
control_titles = ["Thrust","Roll(deg)", "Pitch(deg)", "Yaw(deg)"];
figure
for i = 1:4
    subplot(2, 2, i)
    plot(u_cl(:,i))
    title(control_titles(i))    
end
