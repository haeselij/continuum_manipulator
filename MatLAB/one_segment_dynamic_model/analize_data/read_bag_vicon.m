%% extract data from rosbag vicon
bag = rosbag('bagfiles/vicon_measurement6_unsymmetric_trackenabled_maxpress2lower.bag');

frames = bag.AvailableFrames;
transform_select = select(bag, 'Topic','/segBase/vrpn_client/raw_transform');
transform_time = transform_select.MessageList.Time;
transform_size_vec = size(transform_time);
transform_size = transform_size_vec(1);

pressure_select = select(bag,'Topic','/druecke');
pressure_struct = readMessages(pressure_select,'DataFormat','struct');
pressure_total= cellfun(@(p) double(p.Press), pressure_struct, 'un', 0);
pressure_total = cell2mat(pressure_total');
pressure_time = pressure_select.MessageList.Time; 

pressure_size_vec = size(pressure_time);
pressure_size = pressure_size_vec(1);

length_select = select(bag,'Topic','/laengen');
length_struct = readMessages(length_select,'Dataformat','struct');
length_total = cellfun(@(p) double(p.Length), length_struct, 'un',0);
length_total = cell2mat(length_total');
length_time = length_select.MessageList.Time;
length_size_vec = size(length_time);
length_size = length_size_vec(1);

pressure_time = pressure_time - length_time(1);
length_time = length_time - length_time(1);