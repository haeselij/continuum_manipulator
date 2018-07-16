%% extract data from rosbag without vicon
bag = rosbag('bagfiles/damping_coefficient_identification_2018-07-03-17-37-16.bag');

length_select = select(bag,'Topic','/laengen');
length_struct = readMessages(length_select,'Dataformat','struct');
length_total = cellfun(@(p) double(p.Length), length_struct, 'un',0);
length_total = cell2mat(length_total');
length_time = length_select.MessageList.Time;
length_time = length_time - length_time(1);
length_size_vec = size(length_time);
length_size = length_size_vec(1);
