

%% Plot pressure inputs.
dt = 0.00005;
pressure_time = zeros(1000,1);

for i = 1000
    pressure_time(i) = -0.00005 +dt*i;
end
% 
% figure; 
% plot(pressure_time, pressure_total(2,:), pressure_time, pressure_total(7,:), pressure_time, pressure_total(9,:))
% title('linear increasing pressure input')
% xlabel('Time [s]')
% ylabel('pressure [bar]')
% legend('p1','p2','p3')
% axis([0 0.088 0 0.15])

figure; 
plot(pressure_time, q_(1,:), pressure_time, q_(2,:), pressure_time, q_(3,:))
title('calculated muscle lengths p2 = 0.25')
xlabel('Time [s]')
ylabel('length [m]')
legend('q1','q2','q3')
axis([0 0.5 0 0.6])

% r_b = 0.026;
% 
% for i =1: 1500  
%   
%     theta(i) =  2/3*(sqrt(q_container(1,i)^2 + q_container(2,i)^2 + q_container(3,i)^2 - q_container(1,i)*q_container(2,i) - q_container(1,i)*q_container(3,i) - q_container(2,i)*q_container(3,i)))/r_b;
%     
% end