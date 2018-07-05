function [potential_energy_term] = GetPotentialEnergyTerm(dw_H_1h__dq1, dw_H_1h__dq2, dw_H_1h__dq3, bending_stiffness_derivative)
    syms m g C q_0 A k D_damp q_1_dot q_2_dot q_3_dot v_1_dot v_2_dot v_3_dot phi theta kappa r_b q_1 q_2 q_3 p_1 p_2 p_3 k_spline;

    w_g = [0; 0; -g; 0];     

    potential_energy_q1 = -m*w_g.'*dw_H_1h__dq1(:,4) + (bending_stiffness_derivative(1)); 
    potential_energy_q2 = -m*w_g.'*dw_H_1h__dq2(:,4) + (bending_stiffness_derivative(2)); 
    potential_energy_q3 = -m*w_g.'*dw_H_1h__dq3(:,4) + (bending_stiffness_derivative(3));

    potential_energy_q1 = potential_energy_q1;
    potential_energy_q2 = potential_energy_q2;
    potential_energy_q3 = potential_energy_q3;
 
    potential_energy_term = [potential_energy_q1; potential_energy_q2; potential_energy_q3];
end