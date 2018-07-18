%% Plot each measured and calculated tip position coordinate.
coordinate_size = 14228;
figure;
plot(length_time, r_tip_pos_set(1,:), transform_time(1:coordinate_size), coordinates_measured(1,:))
xlabel('time [s]')
ylabel('length [m]')
set(gca,'fontsize', 12);
legend('x_{calculated}','x_{measured}')
axis([0 150 -0.2 0.15])
%%
figure;
plot(length_time, r_tip_pos_set(2,:), transform_time(1:coordinate_size), coordinates_measured(2,:))
xlabel('time [s]')
ylabel('length [m]')
set(gca,'fontsize', 12);
legend('y_{calculated}','y_{measured}')
axis([0 150 -0.2 0.15])
%%
figure;
plot(length_time, r_tip_pos_set(3,:), transform_time(1:coordinate_size), coordinates_measured(3,1:coordinate_size))
xlabel('time [s]')
ylabel('length [m]')
set(gca,'fontsize', 12);
legend('z_{calculated}','z_{measured}')
axis([0 140 0. 0.35])

%% Plot measured tip position
figure;
plot(transform_time(1:coordinate_size), coordinates_measured(1,:), transform_time(1:coordinate_size), coordinates_measured(2,:), transform_time(1:coordinate_size), coordinates_measured(3,:))
title('measured tip position')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x','y', 'z')
axis([0 150 -0.2 0.4 ])

%% Plot calculated tip position
figure;
plot(length_time, r_tip_pos_set(1,:), length_time,r_tip_pos_set(2,:),length_time,r_tip_pos_set(3,:))
title('calculated tip position')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x','y', 'z')
axis([0 150 -0.2 0.4])


%% Plot angle theta
figure;
plot(length_time, theta(1,:))
title('angle theta')
xlabel('time [s]')
ylabel('theta [rad]')
axis([0 150 0 1.7])

%% Plot angle phi
figure;
plot(length_time, theta(2,:))
title('angle phi')
xlabel('time [s]')
ylabel('phi [rad]')
axis([0 150 -3.2 3.2])

%% Plot pressure inputs.
figure; 
plot(pressure_time, pressure_total(2,:), pressure_time, pressure_total(7,:), pressure_time, pressure_total(9,:))
title('pressure')
xlabel('Time [s]')
ylabel('Length [m]')
legend('p1','p2', 'p3')

%% Plot draw wire lengths
figure;
plot(length_time, length_total(1:3,:))
title('draw wire length')
xlabel('Time [s]')
ylabel('...')
legend('s1', 's2','s3')

%% Plot each measured against each calculated tip position coordinate based on synchronized data
figure;
plot(sync_time_table.Time, sync_time_table.x_r_tip_pos_time_table, sync_time_table.Time, sync_time_table.x_coordinates_measured_time_table)
title('error plot x of synchronized data sets')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x_{calculated}','x_{measured}')
axis([0 145.0 -0.2 0.4])

%%
figure;
plot(sync_time_table.Time, sync_time_table.y_r_tip_pos_time_table, sync_time_table.Time, sync_time_table.y_coordinates_measured_time_table)
title('error plot y of synchronized data sets')
xlabel('Time [s]')
ylabel('Length [m]')
legend('y_{calculated}','y_{measured}')
axis([0 150 -0.2 0.4])

%%
figure;
plot(sync_time_table.Time, sync_time_table.z_r_tip_pos_time_table, sync_time_table.Time, sync_time_table.z_coordinates_measured_time_table)
title('error plot z of synchronized data sets')
xlabel('Time [s]')
ylabel('Length [m]')
legend('z_{calculated}','z_{measured}')
axis([0 145.0 -0.2 0.4])

%% Plot relative norm error
figure;
plot(sync_time_table.Time, absolute_error_norm*1000)
title('absolute norm error plot')
xlabel('Time [s]')
ylabel('Error [mm]')
%axis([0 145 0 15])

%% Plot relative error of each coordinate
figure;
plot(sync_time_table.Time, absolute_error_vector(1,:), sync_time_table.Time, absolute_error_vector(2,:), sync_time_table.Time, absolute_error_vector(3,:))
title('relative error plot of each coordinate')
xlabel('Time [s]')
ylabel('Error [%]')
legend('x_{error}','y_{error}','z_{error}')
axis([0 145.0 0 15])
