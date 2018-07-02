function [q_1, q_2, q_3, v_1, v_2, v_3] = LagrangeEquation(generalized_force_term, kinetic_energy_term, potential_energy_term, q_dot, v_1, v_2, v_3)
syms q_1_n q_2_n q_3_n v_1_n v_2_n v_3_n 

eqn_1 = kinetic_energy_term(1) + potential_energy_term(1) == generalized_force_term(1);
eqn_2 = kinetic_energy_term(2) + potential_energy_term(2) == generalized_force_term(2);
eqn_3 = kinetic_energy_term(3) + potential_energy_term(3) == generalized_force_term(3);

auxilary_eqn_1 = v_1 == q_dot(1);
auxilary_eqn_2 = v_2 == q_dot(2);
auxilary_eqn_3 = v_3 == q_dot(3);

eqns = [eqn_1 eqn_2 eqn_3 auxilary_eqn_1 auxilary_eqn_2 auxilary_eqn_3];
vars = [q_1_n q_2_n q_3_n v_1_n v_2_n v_3_n];
init_guess = [0.3014, 0.3014, 0.3014, 0.25, -0.25, -0.25];
[q_1, q_2, q_3, v_1, v_2, v_3] = vpasolve(eqns, vars, init_guess);

end


