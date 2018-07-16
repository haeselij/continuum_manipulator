% cut data to one damping measurement
j = 1;
m = 0.256;
k = 1108.8;

for i = 1 : length_size
    
    if length_time(i) < 7.15 && length_time(i) > 2.
        
    dw_1(j) = 0.3014 -0.0384/1024*(length_total(1,i)-length_total(2,1));
    dw_2(j) = 0.3014 -0.0384/1024*(length_total(2,i)-length_total(2,1));
    dw_3(j) = 0.3014 -0.0384/1024*(length_total(3,i)-length_total(2,1));

    j = j+1;
    
    end
end

% maximum amplitude points points from measurement
time = [1.0015; 1.377; 1.725; 2.09; 2.438; 2.806; 3.134; 3.507; 3.863;4.211];
time = time-1; %time offset -> set start time to zero

displacement = [0.006812; 0.004562; 0.002875; 0.001862;0.001262; 0.000925; 0.000625; 0.0004375; 0.00025;0];
D_damp = 2*3*1.151*(m*k)^0.5; %values from exponential fitting in matlab (ctftool)
%%
figure;
hold on;
plot(length_time(1:534)-1, dw_1(:) - 0.2962)
plot(time,displacement,'*')
t = 0:0.1:5;
f = 0.0068*exp(-t*1.151);
plot(t,f)
axis([0 5 -0.002 0.007 ])
set(gca,'fontsize', 12);

function [x,y,z,w, theta_1, phi_1] = get_tip(q_11, q_12, q_13, i) 
    if (q_11 == q_12 && q_12 == q_13)
        x = 0;
        y= 0;
        z = .3014;
        w= 1;
        
        theta_1 = 0;
        phi_1 = 0;
        
    else 
        theta_1 = 2/3*(sqrt(q_11^2 + q_12^2 + q_13^2 - q_11*q_12 - q_11*q_13 - q_12*q_13))/0.008;
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
        
        r_tip_pos = R*w_H_1h*[0; 0; 0; 1];
 
        y = - r_tip_pos(1);
        x = r_tip_pos(2);
        z = r_tip_pos(3);
        w = r_tip_pos(4);
        
    end
end