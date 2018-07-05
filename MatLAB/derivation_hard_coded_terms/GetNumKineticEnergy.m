function [kinetic_energy_term_num] = GetNumKineticEnergy(q_1, q_2, q_3, v_1, v_2, v_3, kinetic_energy_term)
syms q_1_n q_2_n q_3_n v_1_n v_2_n v_3_n 

kinetic_energy_term_num = subs(kinetic_energy_term, [q_1_n, q_2_n, q_3_n, v_1_n, v_2_n, v_3_n],[q_1, q_2, q_3, v_1, v_2, v_3]);

end

