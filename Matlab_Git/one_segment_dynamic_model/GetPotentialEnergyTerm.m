function [potential_energy_term] = GetPotentialEnergyTerm(q_1, q_2, q_3, dw_H_1h__dq1, dw_H_1h__dq2, dw_H_1h__dq3, kappa, l_bar, phi, bending_stiffness_derivative)
    w_g = [0; 0; 9.81; 0];
    m = 0.256;
    r_b = 0.035; %[m]
    C_1 = .264; % [Nm/rad]
    C_2 = 2*C_1; %TODO!!!
    C_01 = 0; %TODO!!!
    C_02 = 0; %TODO!!!
    l_bar = (q_1 + q_2 + q_3)/3;
    roundn = @(x,n) round(x*10^n)./10^n;
% if q_1 == q_2 && q_1 == q_3
%     q_1 = q_1 -0.000000001;
%     q_3 = q_3 + 0.000000001;
%     
% end        
%     d_ksys__dq1(1,1) =                                                                                                                                                   -(2*3^(1/2)*(q_2 - q_3))/(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2);
%     d_ksys__dq1(2,1) =  - (2*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_2 - 2*q_1 + q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
%     d_ksys__dq1(3,1) =                                                                                                                                                                                                                 1/3;
% 
% 
% 
%     d_ksys__dq2(1,1) =                                                                         -((3^(1/2)/(q_2 - 2*q_1 + q_3) - (3^(1/2)*(q_2 - q_3))/(q_2 - 2*q_1 + q_3)^2)*(q_2 - 2*q_1 + q_3)^2)/(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2);
%     d_ksys__dq2(2,1) =  - (2*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_1 - 2*q_2 + q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
%     d_ksys__dq2(3,1) =                                                                                                                                                                                                                    1/3;
% 
% 
%     d_ksys__dq3(1,1) =                                                                        ((3^(1/2)/(q_2 - 2*q_1 + q_3) + (3^(1/2)*(q_2 - q_3))/(q_2 - 2*q_1 + q_3)^2)*(q_2 - 2*q_1 + q_3)^2)/(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2);
%     d_ksys__dq3(2,1) = - (2*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_1 + q_2 - 2*q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
%     d_ksys__dq3(3,1) =                                                                                                                                                                                                                   1/3;
% 
% 
% d_M_bend_dq1 = -(C_1*kappa*sin((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2))/3 + (C_2*kappa*cos((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2))/3 - (3*3^(1/2)*C_1*kappa*pi*cos((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2)*(q_2 - q_3)*(q_1/3 + q_2/3 + q_3/3))/(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2) + (3*3^(1/2)*C_2*kappa*pi*sin((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2)*(q_2 - q_3)*(q_1/3 + q_2/3 + q_3/3))/(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2)
% d_M_bend_dq2 = -(C_1*kappa*sin((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2))/3 + (C_2*kappa*cos((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2))/3 - (3*C_1*kappa*pi*cos((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2)*(3^(1/2)/(q_2 - 2*q_1 + q_3) - (3^(1/2)*(q_2 - q_3))/(q_2 - 2*q_1 + q_3)^2)*(q_2 - 2*q_1 + q_3)^2*(q_1/3 + q_2/3 + q_3/3))/(2*(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2)) + (3*C_2*kappa*pi*sin((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2)*(3^(1/2)/(q_2 - 2*q_1 + q_3) - (3^(1/2)*(q_2 - q_3))/(q_2 - 2*q_1 + q_3)^2)*(q_2 - 2*q_1 + q_3)^2*(q_1/3 + q_2/3 + q_3/3))/(2*(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2))
% d_M_bend_dq3 = -(C_1*kappa*sin((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2))/3 + (C_2*kappa*cos((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2))/3 + (3*C_1*kappa*pi*cos((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2)*(3^(1/2)/(q_2 - 2*q_1 + q_3) + (3^(1/2)*(q_2 - q_3))/(q_2 - 2*q_1 + q_3)^2)*(q_2 - 2*q_1 + q_3)^2*(q_1/3 + q_2/3 + q_3/3))/(2*(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2)) - (3*C_2*kappa*pi*sin((3*pi*atan2(-3^(1/2)*(q_2 - q_3), q_2 - 2*q_1 + q_3))/2)*(3^(1/2)/(q_2 - 2*q_1 + q_3) + (3^(1/2)*(q_2 - q_3))/(q_2 - 2*q_1 + q_3)^2)*(q_2 - 2*q_1 + q_3)^2*(q_1/3 + q_2/3 + q_3/3))/(2*(3*(q_2 - q_3)^2 + (q_2 - 2*q_1 + q_3)^2))
%  
% if phi == 0
%     phi = 0.001;
% end
%     d_M_bend_dksys = [ (C_02 + C_2*kappa*l_bar)*(cos(phi)/(sin(phi)^2*(1/sin(phi) - 1/sin(phi - pi/3))) - (cos(phi)/sin(phi)^2 - cos(phi - pi/3)/sin(phi - pi/3)^2)/(sin(phi)*(1/sin(phi) - 1/sin(phi - pi/3))^2)) - (cos(phi)*(C_01 + C_1*kappa*l_bar))/(sin(phi)^2*(1/sin(phi) - 1/sin(phi - pi/3))) + ((C_01 + C_1*kappa*l_bar)*(cos(phi)/sin(phi)^2 - cos(phi - pi/3)/sin(phi - pi/3)^2))/(sin(phi)*(1/sin(phi) - 1/sin(phi - pi/3))^2), (C_1*l_bar)/(sin(phi)*(1/sin(phi) - 1/sin(phi - pi/3))) - C_2*l_bar*(1/(sin(phi)*(1/sin(phi) - 1/sin(phi - pi/3))) - 1), (C_1*kappa)/(sin(phi)*(1/sin(phi) - 1/sin(phi - pi/3))) - C_2*kappa*(1/(sin(phi)*(1/sin(phi) - 1/sin(phi - pi/3))) - 1)];
%   
if q_1 == q_2 && q_2 ==  q_3
    potential_energy_q1 = 0;
    potential_energy_q2 = 0;
    potential_energy_q3 = 0;
else
    potential_energy_q1 = m*w_g.'*dw_H_1h__dq1(:,4) ;%+ bending_stiffness_derivative(1); %C*l_bar%% WIE IST HIER DAS VORZEICHEN??
    potential_energy_q2 = m*w_g.'*dw_H_1h__dq2(:,4) ;%+ bending_stiffness_derivative(2); %+ d_M_bend_dq1;
    potential_energy_q3 = m*w_g.'*dw_H_1h__dq3(:,4) ;%+ bending_stiffness_derivative(3);% + d_M_bend_dq3;

    potential_energy_q1 = roundn(potential_energy_q1,10);
    potential_energy_q2 = roundn(potential_energy_q2,10);
    potential_energy_q3 = roundn(potential_energy_q3,10);
 end
    potential_energy_term = [potential_energy_q1; potential_energy_q2; potential_energy_q3];
end