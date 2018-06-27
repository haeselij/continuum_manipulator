function [q_1, q_2 q_3] = LagrangeEquation(dt, generalized_force_term,kinetic_energy_term, potential_energy_term)
syms q_1 q_2 q_3 v_1 v_2 v_3 
syms d_q_1 d_q_2 d_q_3 d2_q_1 d2_q_2 d2_q_3

kinetic_energy_t = simplify(kinetic_energy_term);
potential_energy_t = simplify(potential_energy_term);
generalized_force_t = simplyfiy(generalized_force_term);
eqn_1 = kinetic_energy_t(1)+ potential_energy_t(1) == generalized_force_t(1);

eqns = [eqn_1 eqn_12 eqn_13 eqn_21 eqn_22 eqn_23 h_eqn_11 h_eqn_12 h_eqn_13 h_eqn_21 h_eqn_22 h_eqn_23];
vars = [x_11_2(t); x_12_2(t); x_13_2(t); x_21_2(t); x_22_2(t); x_23_2(t); q_11(t); q_12(t); q_13(t); q_21(t); q_22(t); q_23(t)];

end


