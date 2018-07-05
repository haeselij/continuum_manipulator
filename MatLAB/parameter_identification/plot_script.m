%% Plot each measured and calculated tip position coordinate.
figure;
plot(length_time, r_tip_pos_set(1,:), transform_time(1:coordinate_size), coordinates_measured(1,:))
title('error plot x')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x_{calculated}','x_{measured}')
axis([0 500 -0.2 0.4])

figure;
plot(length_time, r_tip_pos_set(2,:), transform_time(1:coordinate_size), coordinates_measured(2,:))
title('error plot y')
xlabel('Time [s]')
ylabel('Length [m]')
legend('y_{calculated}','y_{measured}')
axis([0 500 -0.2 0.4])

figure;
plot(length_time, r_tip_pos_set(3,:), transform_time(1:coordinate_size), coordinates_measured(3,:))
title('error plot z')
xlabel('Time [s]')
ylabel('Length [m]')
legend('z_{calculated}','z_{measured}')
axis([0 500 -0.2 0.4])

%% Plot measured tip position
figure;
plot(transform_time(1:coordinate_size), coordinates_measured(1,:), transform_time(1:coordinate_size), coordinates_measured(2,:), transform_time(1:coordinate_size), coordinates_measured(3,:))
title('measured tip position')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x','y', 'z')
axis([0 300 -0.2 0.4 ])

%% Plot calculated tip position
figure;
plot(length_time, r_tip_pos_set(1,:), length_time,r_tip_pos_set(2,:),length_time,r_tip_pos_set(3,:))
title('calculated tip position')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x','y', 'z')
axis([0 300 -0.2 0.4])


%% Plot angle theta
figure;
plot(length_time, theta(1,:))
title('angle theta')
xlabel('time [s]')
ylabel('theta [rad]')
axis([0 4500 0 1.7])

%% Plot angle phi
figure;
plot(length_time, theta(2,:))
title('angle phi')
xlabel('time [s]')
ylabel('phi [rad]')
axis([0 4500 -3.2 3.2])

%% Plot pressure inputs.
figure; 
plot(pressure_time, pressure_total(1,:), pressure_time, pressure_total(2,:), pressure_time, pressure_total(3,:))
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
axis([0 200 200 600])
%% visualization of the results

% figure;
% 
% scatter3( x_measured, y_measured, z_measured ,40)
% title('measured tip position')
% cd =colorbar;
% cd.Label.String ='Time [s]';
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')

%% Plot each measured against each calculated tip position coordinate based on synchronized data
figure;
plot(sync_time_table.Time, sync_time_table.x_r_tip_pos_time_table, sync_time_table.Time, sync_time_table.x_coordinates_measured_time_table)
title('error plot x of synchronized data sets')
xlabel('Time [s]')
ylabel('Length [m]')
legend('x_{calculated}','x_{measured}')
axis([0 4500 -0.2 0.4])

%%
figure;
plot(sync_time_table.Time, sync_time_table.y_r_tip_pos_time_table, sync_time_table.Time, sync_time_table.y_coordinates_measured_time_table)
title('error plot y of synchronized data sets')
xlabel('Time [s]')
ylabel('Length [m]')
legend('y_{calculated}','y_{measured}')
axis([0 4500 -0.2 0.4])

%%
figure;
plot(sync_time_table.Time, sync_time_table.z_r_tip_pos_time_table, sync_time_table.Time, sync_time_table.z_coordinates_measured_time_table)
title('error plot z of synchronized data sets')
xlabel('Time [s]')
ylabel('Length [m]')
legend('z_{calculated}','z_{measured}')
axis([0 4500 -0.2 0.4])

%% Plot relative norm error
figure;
plot(sync_time_table.Time, relative_norm_error)
title('relative norm error plot')
xlabel('Time [s]')
ylabel('Error [%]')
axis([0 4500 0 15])

%% Plot relative error of each coordinate
figure;
plot(sync_time_table.Time, relative_x_coordinate_error, sync_time_table.Time, relative_y_coordinate_error, sync_time_table.Time, relative_z_coordinate_error)
title('relative error plot of each coordinate')
xlabel('Time [s]')
ylabel('Error [%]')
legend('x_{error}','y_{error}','z_{error}')
axis([0 4500 0 15])

%% Plot interpolated draw wire sensor data
figure;
plot(interpolated_time_table.pressure_time_converted, interpolated_time_table.s_1, interpolated_time_table.pressure_time_converted, interpolated_time_table.s_2, interpolated_time_table.pressure_time_converted, interpolated_time_table.s_3)
title('interpolated draw wire data')
xlabel('Time [s]')
ylabel('... [%]')
legend('s_{1}','s_{2}','s_{3}')

%% Plot normalized pressure and draw wire sensor
figure;
plot(interpolated_time_table_fast.pressure_time_converted, length_normalized(1,:), interpolated_time_table_fast.pressure_time_converted, length_normalized(2,:), interpolated_time_table_fast.pressure_time_converted, length_normalized(3,:), interpolated_time_table_fast.pressure_time_converted, pressure_normalized)
title('interpolated pressure data')
xlabel('Time [s]')
ylabel('Pressure [bar]')
legend('p_1','p_2','p_3')