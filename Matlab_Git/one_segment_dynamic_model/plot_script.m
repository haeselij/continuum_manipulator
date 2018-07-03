

%% Plot pressure inputs.
dt = 0.001;
pressure_time = zeros(100,i);

for i = 1:88
    pressure_time(i+1) = pressure_time(i) + dt;
end

figure; 
plot(pressure_time, pressure_total(2,:), pressure_time, pressure_total(7,:), pressure_time, pressure_total(9,:))
title('linear increasing pressure input')
xlabel('Time [s]')
ylabel('pressure [bar]')
legend('p1','p2','p3')
axis([0 0.088 0 0.15])

figure; 
plot(pressure_time, q_container(1,:), pressure_time, q_container(2,:), pressure_time, q_container(3,:))
title('calculated muscle lengths')
xlabel('Time [s]')
ylabel('length [m]')
legend('q1','q2','q3')
axis([0 0.088 0 0.4])


