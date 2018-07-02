function [kinetic_energy_term] = GetKineticEnergyTerm(dw_H_1h__dq1,dw_H_1h__dq2,dw_H_1h__dq3, w_H_2dot)
    m = 0.256;

    kinetic_energy_term_1 = w_H_2dot.'*dw_H_1h__dq1(:,4);
    kinetic_energy_term_2 = w_H_2dot.'*dw_H_1h__dq2(:,4);
    kinetic_energy_term_3 = w_H_2dot.'*dw_H_1h__dq3(:,4);
    
    kinetic_energy_term = m*[kinetic_energy_term_1; kinetic_energy_term_2; kinetic_energy_term_3];
end