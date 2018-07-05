function [eqns] = LagrangeEquation(generalized_force_term, kinetic_energy_term, potential_energy_term)
    eqn_1 = kinetic_energy_term(1) + potential_energy_term(1) == generalized_force_term(1);
    eqn_2 = kinetic_energy_term(2) + potential_energy_term(2) == generalized_force_term(2);
    eqn_3 = kinetic_energy_term(3) + potential_energy_term(3) == generalized_force_term(3);

    eqns = [eqn_1 eqn_2 eqn_3];
end


