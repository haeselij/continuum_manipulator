function [potential_energy_term] = GetPotentialEnergyTerm()
   dkappa__dq1 = - (2*(q_1^2 - q_1*q_3 - q_2*q_3 - q_1*q_2 + q_2^2 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_2 - 2*q_1 + q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
   dkappa__dq2 = - (2*(q_1^2 - q_1*q_3 - q_2*q_3 - q_1*q_2 + q_2^2 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_1 - 2*q_2 + q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
   dkappa__dq3 = - (2*(q_1^2 - q_1*q_3 - q_2*q_3 - q_1*q_2 + q_2^2 + q_3^2)^(1/2))/(9*r_b*(q_1/3 + q_2/3 + q_3/3)^2) - (q_1 + q_2 - 2*q_3)/(3*r_b*(q_1/3 + q_2/3 + q_3/3)*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
   
   potential_energy_q1 = m*w_g.'*dw_H_1h__dq1(:,4) + C*l_bar*dkappa__dq1;
   potential_energy_q2 =
   potential_energy_q3
end