%% Visualization of Kinematic Model
%start roscore
%run continuum_manipulator_visualization
%open rviz

r_b = 0.026;

q_11_set = length_total(1,:);
q_12_set = length_total(2,:);
q_13_set = length_total(3,:);



for i = 1 : length_size
    
    dw_1 = 0.3014 +0.0384/1024*(q_11_set(i)-q_11_set(1));
    dw_2 = 0.3014 +0.0384/1024*(q_12_set(i)-q_12_set(1));
    dw_3 = 0.3014 +0.0384/1024*(q_13_set(i)- q_13_set(1));
    [x_, theta_, phi_] =  GetXThetaPhi(dw_1, dw_2, dw_3);
    [q_(1,i), q_(2,i), q_(3,i)] = GetBalgLength(theta_, phi_);
 
end


% ros publisher 

rosinit 
node = robotics.ros.Node('/muscle');
rate = robotics.ros.Rate(node,240);
pub = rospublisher('/lengths','std_msgs/Float32MultiArray');
msg = rosmessage(pub);

reset(rate);

for i= 1:length_size
    
   msg.Data(1) =   q_(1,i);
   msg.Data(2) =   q_(2,i);
   msg.Data(3) =   q_(3,i);
   
   send(pub,msg);
   waitfor(rate);
   
end
rosshutdown
%%

function [x,theta, phi] = GetXThetaPhi(q_1, q_2, q_3)
    theta = (2/3)*((sqrt(q_1^2 + q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3 - q_2*q_3))/0.01);
    phi = atan2((sqrt(3)*(q_3 - q_2)),(q_2 + q_3 - 2*q_1));
    l_bar = (q_1 + q_2 + q_3)/3.0;
    kappa = theta/l_bar;
    
    x = cos(phi)*(1-cos(theta))/kappa;
end

function [q_11, q_12, q_13] = GetBalgLength(theta, phi)
    phi_0 = pi/3;
    l_bar_1 = 0.3014;
    r_b = 0.026;
    kappa = theta/l_bar_1;
    if theta == 0
        q_11 = 0.3014;
        q_12 = 0.3014;
        q_13 = 0.3014;
    else
    q_11 = theta*(1/kappa - r_b*cos(phi + phi_0));
    q_12 = theta*(1/kappa - r_b*cos(phi + 2*pi/3 + phi_0));
    q_13 = theta*(1/kappa - r_b*cos(phi + 4*pi/3 + phi_0));
    end
end
