

%% Plot pressure inputs.
dt = 0.0005;
pressure_time = zeros(640,1);

for i = 1:640
    pressure_time(i) = -0.01 +dt*i;
end

% figure; 
% plot(pressure_time, pressure_total(2,:), pressure_time, pressure_total(7,:), pressure_time, pressure_total(9,:))
% title('linear increasing pressure input')
% xlabel('Time [s]')
% ylabel('pressure [bar]')
% legend('p1','p2','p3')
% axis([0 0.088 0 0.15])

figure; 
plot(pressure_time, q_container(1,1:640), pressure_time, q_container(2,1:640), pressure_time, q_container(3,1:640))
title('calculated muscle lengths')
xlabel('Time [s]')
ylabel('length [m]')
legend('q1','q2','q3')
%axis([0 0.4 0 0.4])


