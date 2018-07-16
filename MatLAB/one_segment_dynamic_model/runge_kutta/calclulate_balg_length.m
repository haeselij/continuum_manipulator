%% Visualization of Kinematic Model
% syms q_1 q_2 q_3
r_b = 0.026;
q_11_set = length_total(1,:);
q_12_set = length_total(2,:);
q_13_set = length_total(3,:);
% w_H_t_total = zeros(4,4,length_size);


for i = 1 : length_size
    dw_1 = 0.3014 +0.0384/1024*(q_11_set(i)-q_11_set(1));
    dw_2 = 0.3014 +0.0384/1024*(q_12_set(i)-q_12_set(1));
    dw_3 = 0.3014 +0.0384/1024*(q_13_set(i)- q_13_set(1));
    l_bar(i) = (dw_1 + dw_2 +dw_3)/3;
 [x_(i), theta_(i), phi_(i)] =  GetXThetaPhi(dw_1, dw_2, dw_3);
 [q_exp(1,i), q_exp(2,i), q_exp(3,i)] = GetBalgLength(theta_(i), phi_(i), l_bar(i));
 
 
 end
%%
% %length_time = length_time -1;
 length_time_1 = length_time -8.07;
%length_time_2 = length_time -32.82;
 length_time_3 = length_time -78.29;
% 
% %%
figure;
plot(length_time_1,q_exp(2,:),length_time_1, q_exp(1,:),length_time_1, q_exp(3,:))
%title('experimental muscle lengths')
xlabel('time [s]')
ylabel('length [m]')
legend('q_1','q_2','q_3')
set(gca,'fontsize', 12);
axis([0 25 0.25 0.35])







%%
function [x,theta, phi] = GetXThetaPhi(q_1, q_2, q_3)
    theta = (2/3)*((sqrt(q_1^2 + q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3 - q_2*q_3))/0.01);
    phi = atan2((sqrt(3)*(q_3 - q_2)),(q_2 + q_3 - 2*q_1));
    l_bar = (q_1 + q_2 + q_3)/3.0;
    kappa = theta/l_bar;
    
    x = cos(phi)*(1-cos(theta))/kappa;
end

function [q_11, q_12, q_13] = GetBalgLength(theta, phi, l_bar_1)
    phi_0 = pi/3;
   % l_bar_1 = 0.3014;
    r_b = 0.026;
    kappa = theta/l_bar_1;
    if theta == 0
        q_11 = 0.3014;
        q_12 = 0.3014;
        q_13 = 0.3014;
    else
    q_11 = theta*(1/kappa - r_b*cos(phi + phi_0));
    q_12 = theta*(1/kappa - r_b*cos(phi + 4*pi/3 + phi_0));
    q_13 = theta*(1/kappa - r_b*cos(phi + 2*pi/3 + phi_0));
    end
end
% %%
% rosinit 
% node = robotics.ros.Node('/muscle');
% rate = robotics.ros.Rate(node,240);
% pub = rospublisher('/lengths','std_msgs/Float32MultiArray');
% msg = rosmessage(pub);
% 
% reset(rate);
% for i= 1:length_size
%    msg.Data(1) =   q_(1,i);
%    msg.Data(2) =   q_(2,i);
%    msg.Data(3) =   q_(3,i);
%    send(pub,msg);
%    waitfor(rate);
% end
% rosshutdown