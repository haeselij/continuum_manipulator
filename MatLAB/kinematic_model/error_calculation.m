%% Derive all time differences
for i = 0:length_size/2-1
    delta_time_lenght(i+1) = length_time(2*i + 2) - length_time(2*i + 1);
end

for i = 0:transform_size/2-1
    delta_time_tf(i+1) = transform_time(2*i + 2) - transform_time(2*i + 1);
end

%% Cut transform_time to the last time of length_time
check_variable = length_time(length_size) < transform_time(transform_size);
cut_position = 0;

if(check_variable)
    fprintf('Last transform_time is bigger. \n')
    
    for i = 1:transform_size -2 
        if(length_time(length_size) > transform_time(i))
            cut_position = i;
            transform_time_cut(i) = transform_time(i);
            coordinates_measured_cut(:,i) = coordinates_measured(:,i);
        end
    end
    transform_frequency = cut_position/transform_time_cut(cut_position);
    length_frequency = length_size/length_time(length_size);
    
    length_time_cut = length_time;
    r_tip_pos_set_cut = r_tip_pos_set;
       
    
else
    fprintf('Last transform_time is bigger. \n')
    for i = 1:length_size
        if(transform_time(transform_size) > length_time(i))
            cut_position = i;
            length_time_cut(i) = length_time(i);
            r_tip_pos_set_cut(:,i) = r_tip_pos_set(:,i);
        end
    length_frequency = cut_position/length_time_cut(cut_position);
    transform_frequency = transform_size/transform_time(transform_size);
    
    transform_time_cut = transform_time;
    coordinates_measured_cut = coordinates_measured;
      
    end
    
end
 
%% Synchronize and interpolate data
coordinates_measured_cut(:,coordinate_size+1) = coordinates_measured_cut(:,coordinate_size);
coordinates_measured_cut(:,coordinate_size+2) = coordinates_measured_cut(:,coordinate_size);

length_time_converted = seconds(length_time_cut.');
transform_time_converted = seconds(transform_time_cut.');
r_tip_pos_time_table = timetable(length_time_converted, r_tip_pos_set_cut(1,:).', r_tip_pos_set_cut(2,:).', r_tip_pos_set_cut(3,:).');
r_tip_pos_time_table.Properties.VariableNames = {'x','y','z'};
coordinates_measured_time_table = timetable(transform_time_converted(1:14229).', coordinates_measured_cut(1,:).', coordinates_measured_cut(2,:).', coordinates_measured_cut(3,:).'); 
coordinates_measured_time_table.Properties.VariableNames = {'x','y','z'};

%sync_time_table = synchronize(coordinates_measured_time_table, r_tip_pos_time_table, 'union', 'linear');
sync_time_table = synchronize(coordinates_measured_time_table, r_tip_pos_time_table, 'regular', 'linear', 'SamplingRate', 100);

%% Calculate error
absolute_error_vector(1,:) = sync_time_table.x_coordinates_measured_time_table - sync_time_table.x_r_tip_pos_time_table;
absolute_error_vector(2,:) = sync_time_table.y_coordinates_measured_time_table - sync_time_table.y_r_tip_pos_time_table;
absolute_error_vector(3,:) = sync_time_table.z_coordinates_measured_time_table - sync_time_table.z_r_tip_pos_time_table;


for i = 1:size(absolute_error_vector,2)
    absolute_error_norm(i) = norm(absolute_error_vector(:,i));
    coordinates_measured_norm(i) = norm([sync_time_table.x_coordinates_measured_time_table(i); sync_time_table.y_coordinates_measured_time_table(i); sync_time_table.z_coordinates_measured_time_table(i)]);
    
    relative_norm_error(i) = 100*absolute_error_norm(i)/coordinates_measured_norm(i);
    
    relative_x_coordinate_error(i) = norm(100*absolute_error_vector(1,i)/sync_time_table.x_coordinates_measured_time_table(i));
    relative_y_coordinate_error(i) = norm(100*absolute_error_vector(2,i)/sync_time_table.y_coordinates_measured_time_table(i));
    relative_z_coordinate_error(i) = norm(100*absolute_error_vector(3,i)/sync_time_table.z_coordinates_measured_time_table(i));
    
end
absolute_error_mean = mean(absolute_error_norm)
std = std(absolute_error_norm)
%% Calculate theta and phi from measured tip position


