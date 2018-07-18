
%% main subroutine, tip positions are calculated here
q_11_set = length_total(1,:);
q_12_set = length_total(2,:);
q_13_set = length_total(3,:);
w_H_t_total = zeros(4,4,length_size);

for i = 1 : length_size
    a = 0.3014 + 0.0384/1024*(q_11_set(i)-q_11_set(1));
    b = 0.3014 + 0.0384/1024*(q_12_set(i)-q_12_set(1));
    c = 0.3014 + 0.0384/1024*(q_13_set(i)- q_13_set(1));
    [x, y, z, w, theta_1, phi_1] = get_tip(a,b,c,i);
    r_tip_pos_set(1,i) = -x;
    r_tip_pos_set(2,i) = -y;
    r_tip_pos_set(3,i) = z;
    theta(1,i)= theta_1;
    theta(2,i) = phi_1; 
end

%% define function
function [x,y,z,w, theta_1, phi_1] = get_tip(q_11, q_12, q_13, i) 
    if (q_11 == q_12 && q_12 == q_13)
        x = 0;
        y= 0;
        z = .3014;
        w= 1;
        
        theta_1 = 0;
        phi_1 = 0;
        
    else 
        theta_1 = 2/3*(sqrt(q_11^2 + q_12^2 + q_13^2 - q_11*q_12 - q_11*q_13 - q_12*q_13))/0.0082;
        l_bar_1 = (q_11 + q_12 + q_13)/3;
        phi_1 = atan2((sqrt(3)*(q_13 - q_12)),(q_12 + q_13 - 2*q_11));
        kappa_1 = theta_1/l_bar_1;

        w_H_1h_1 = [cos(phi_1)^2*(cos(theta_1) - 1) + 1, sin(phi_1)*cos(phi_1)*(cos(theta_1) - 1), cos(phi_1)*sin(theta_1), (cos(phi_1)*(1 - cos(theta_1)))/kappa_1];
        w_H_1h_2 = [sin(phi_1)*cos(phi_1)*(cos(theta_1) - 1), cos(phi_1)^2*(1 - cos(theta_1)) + cos(theta_1), sin(phi_1)*sin(theta_1), (sin(phi_1)*(1 - cos(theta_1)))/kappa_1];
        w_H_1h_3 = [-cos(phi_1)*sin(theta_1), -sin(phi_1)*sin(theta_1), cos(theta_1), sin(theta_1)/kappa_1];
        w_H_1h_4 = [0, 0, 0, 1];
        w_H_1h = [w_H_1h_1; w_H_1h_2; w_H_1h_3; w_H_1h_4];
        
        
        alpha = pi*-17.5/180;
        R = [cos(alpha), -sin(alpha), 0, 0; sin(alpha), cos(alpha), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
        
        w_H_t_total(:,:,i) = R*w_H_1h;
        r_tip_pos = R*w_H_1h*[0; 0; 0; 1];
 
        y = - r_tip_pos(1);
        x = r_tip_pos(2);
        z = r_tip_pos(3);
        w = r_tip_pos(4);
    end
end