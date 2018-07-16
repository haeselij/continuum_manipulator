%% interpolate data
length_time_converted = seconds(length_time);
pressure_time_converted = seconds(pressure_time);

length_time_table = timetable(length_time_converted, length_total(1,:).', length_total(2,:).', length_total(3,:).');
length_time_table.Properties.VariableNames = {'s_1', 's_2', 's_3'};

pressure_time_table = timetable(pressure_time_converted, pressure_total(1,:).', pressure_total(2,:).', pressure_total(3,:).');
pressure_time_table.Properties.VariableNames = {'p1', 'p2', 'p3'};


interpolated_time_table_fast = synchronize(pressure_time_table, length_time_table, 'regular', 'linear', 'SamplingRate', 200);

%% Filter static points
num_data_points_to_even = 1062; % experience value read out from bags
global_iterator = 1;
num_of_presssure_increments = 9;
all_static_positions_stored = false;

 for j = 1:size(interpolated_time_table_fast,1)
     
        if(interpolated_time_table_fast.p1(j) > 0.1 || interpolated_time_table_fast.p2(j) > 0.1 | interpolated_time_table_fast.p3(j) > 0.1)
            
            for i = 0:num_of_presssure_increments
                
                starting_position = pos_of_first_pressure + i*2000 + num_data_points_to_even;
                end_position = pos_of_first_pressure + (i + 1)*2000;
                new_global_iterator = global_iterator + end_position - starting_position;
                
                length_static(1, global_iterator:new_global_iterator) = interpolated_time_table_fast.s_1(starting_position:end_position).'; 
                length_static(2, global_iterator:new_global_iterator) = interpolated_time_table_fast.s_2(starting_position:end_position).'; 
                length_static(3, global_iterator:new_global_iterator) =  interpolated_time_table_fast.s_3(starting_position:end_position).';
                pressure_static(1, global_iterator:new_global_iterator) = interpolated_time_table_fast.p1(starting_position:end_position).';
                pressure_static(2, global_iterator:new_global_iterator) = interpolated_time_table_fast.p2(starting_position:end_position).'; 
                pressure_static(3, global_iterator:new_global_iterator) = interpolated_time_table_fast.p3(starting_position:end_position).';
            
                global_iterator = new_global_iterator;
                
            end
            
            all_static_positions_stored = true;
            
        else
            
            pos_of_first_pressure = j;
        
        end
        
        if(all_static_positions_stored)
            
            break
            
        end
 end
 
%% main routine
initial_length = [interpolated_time_table_fast.s_1(1); interpolated_time_table_fast.s_2(1); interpolated_time_table_fast.s_3(1)];

for i = 1:size(length_static,2)
    
        length_effective(:, i) = 0.3014 + (0.0384/1024)*(length_static(:, i) - initial_length);
    
        [q_(:, i), theta_draw_wire_(i, 1), phi_muslce_(i, 1), kappa_draw_wire_(i, 1)] = GetBalgLength(length_effective(:,i));
        
end

%% Define needed functions
function [q, theta, phi_muscle, kappa] = GetBalgLength(s)
    phi_0 = pi/3;
    r_b = 0.026; 
    r_s = 0.010;
    
    theta = (2/3)*((sqrt(s(1)^2 + s(2)^2 + s(3)^2 - s(1)*s(2) - s(1)*s(3) - s(2)*s(3)))/r_s);
    phi = atan2((sqrt(3)*(s(3) - s(2))), (s(2) + s(3) - 2*s(1)));
    phi_muscle = phi + pi/3;
    l_bar = (s(1) + s(2) + s(3))/3;
    kappa = theta/l_bar;
   
    if theta == 0
        
        q(1, 1) = 0.3014;
        q(2, 1) = 0.3014;
        q(3, 1) = 0.3014;
        
    else
        
        q(1, 1) = theta*(1/kappa - r_b*cos(phi + phi_0));
        q(2, 1) = theta*(1/kappa - r_b*cos(phi + 4*pi/3 + phi_0));
        q(3, 1) = theta*(1/kappa - r_b*cos(phi + 2*pi/3 + phi_0));
        
    end
end


