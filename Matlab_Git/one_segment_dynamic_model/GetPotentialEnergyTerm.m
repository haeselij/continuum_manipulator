function [potential_energy_term] = GetPotentialEnergyTerm(q_1, q_2, q_3, dw_H_1h__dq1,dw_H_1h__dq2, dw_H_1h__dq3)
    w_g = [0; 0; 9.81;0];
    m = 0.256;
    r_b = 0.035; %[m]
    C = 2.4; % [Nm/rad]
    l_bar = (q_1 + q_2 + q_3)/3;

    dkappa__dq1 = - (2*(q_1^2 - q_1*q_3 - q_2*q_3 - q_1*q_2 + q_2^2 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_2 - 2*q_1 + q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
    dkappa__dq2 = - (2*(q_1^2 - q_1*q_3 - q_2*q_3 - q_1*q_2 + q_2^2 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_1 - 2*q_2 + q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
    dkappa__dq3 = - (2*(q_1^2 - q_1*q_3 - q_2*q_3 - q_1*q_2 + q_2^2 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_1 + q_2 - 2*q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));

    potential_energy_q1 = m*w_g.'*dw_H_1h__dq1(:,4) + C*l_bar*dkappa__dq1;
    potential_energy_q2 = m*w_g.'*dw_H_1h__dq2(:,4) + C*l_bar*dkappa__dq2;
    potential_energy_q3 =  m*w_g.'*dw_H_1h__dq3(:,4) + C*l_bar*dkappa__dq3;

    potential_energy_term = [potential_energy_q1; potential_energy_q2; potential_energy_q3];
end