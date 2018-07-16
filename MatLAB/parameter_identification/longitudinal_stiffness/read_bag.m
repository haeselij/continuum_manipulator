%% extract data from rosbag without vicon
bag = rosbag('bagfiles/stiffness_measurement_2_2018-06-26-14-37-43.bag');

length_select = select(bag,'Topic','/laengen');
length_struct = readMessages(length_select,'Dataformat','struct');
length_total = cellfun(@(p) double(p.Length), length_struct, 'un',0);
length_total = cell2mat(length_total');
length_time = length_select.MessageList.Time;
length_size_vec = size(length_time);
length_size = length_size_vec(1);
time_zero = length_time(1);
length_time = length_time - length_time(1);


pressure_select = select(bag,'Topic','/druecke');
pressure_struct = readMessages(pressure_select,'DataFormat','struct');
pressure_total= cellfun(@(p) double(p.Press), pressure_struct, 'un', 0);
pressure_total = cell2mat(pressure_total');
pressure_time = pressure_select.MessageList.Time; 
pressure_time = pressure_time - time_zero;
pressure_size_vec = size(pressure_time);
pressure_size = pressure_size_vec(1);