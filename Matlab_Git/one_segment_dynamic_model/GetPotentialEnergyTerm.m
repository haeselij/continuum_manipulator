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

    potential_energy_q1 = m*w_g.'*dw_H_1h__dq1(:,4) ;%+ real(bending_stiffness_derivative(1)); %C*l_bar%% WIE IST HIER DAS VORZEICHEN??
    potential_energy_q2 = m*w_g.'*dw_H_1h__dq2(:,4) ;%+ real(bending_stiffness_derivative(2)); %+ d_M_bend_dq1;
    potential_energy_q3 = m*w_g.'*dw_H_1h__dq3(:,4) ;%+ real(bending_stiffness_derivative(3));% + d_M_bend_dq3;

    potential_energy_q1 = roundn(potential_energy_q1,10);
    potential_energy_q2 = roundn(potential_energy_q2,10);
    potential_energy_q3 = roundn(potential_energy_q3,10);
 
    potential_energy_term = [potential_energy_q1; potential_energy_q2; potential_energy_q3];
end