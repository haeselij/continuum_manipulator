function [kinetic_energy_term] = GetKineticEnergyTerm(dw_H_1h__dq1, dw_H_1h__dq2, dw_H_1h__dq3, w_H_2dot)

    syms m C q_0 A k D_damp q_1_dot q_2_dot q_3_dot v_1_dot v_2_dot v_3_dot phi theta kappa r_b q_1 q_2 q_3 p_1 p_2 p_3 k_spline;

    kinetic_energy_term_1 = w_H_2dot.'*dw_H_1h__dq1(:,4);
    kinetic_energy_term_2 = w_H_2dot.'*dw_H_1h__dq2(:,4);
    kinetic_energy_term_3 = w_H_2dot.'*dw_H_1h__dq3(:,4);

    kinetic_energy_term = m*[kinetic_energy_term_1; kinetic_energy_term_2; kinetic_energy_term_3];

end
