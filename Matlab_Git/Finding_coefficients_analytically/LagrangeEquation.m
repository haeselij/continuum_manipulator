function [q_1, q_2, q_3, v_1, v_2, v_3] = LagrangeEquation(generalized_force_term, kinetic_energy_term, potential_energy_term, q_dot, v_1, v_2, v_3, q_1, q_2, q_3)
syms C_1 D_damp k_spline

eqn_1 = kinetic_energy_term(1) + potential_energy_term(1) == generalized_force_term(1);
eqn_2 = kinetic_energy_term(2) + potential_energy_term(2) == generalized_force_term(2);
eqn_3 = kinetic_energy_term(3) + potential_energy_term(3) == generalized_force_term(3);

auxilary_eqn_1 = v_1 == q_dot(1);
auxilary_eqn_2 = v_2 == q_dot(2);
auxilary_eqn_3 = v_3 == q_dot(3);

vars = [k_spline C_1 D_damp];
eqns = [eqn_1 eqn_2 eqn_3, auxilary_eqn_1, auxilary_eqn_2,auxilary_eqn_3];
eqns = simplify(eqns);

init_guess = [1000, 0.3, 10];

[k_spline C_1 D_damp] = solve(eqns, vars);
end


