
   [w_H_1h_, w_H_2dot_, dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_]= SetMatrices();
  
   [generalized_force_term_] = GetGeneralizedForceTerm(w_H_1h_, dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_);
   
   [bending_stiffness_derivative_] = GetBendingStiffnessDeriv();
   [kinetic_energy_term_] = GetKineticEnergyTerm(dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, w_H_2dot_);
   [potential_energy_term_] = GetPotentialEnergyTerm(dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, bending_stiffness_derivative_);
 
   [eqns_] = LagrangeEquation(generalized_force_term_, kinetic_energy_term_, potential_energy_term_);

   %%
    str1 = sprintf('eqn1_singularität');
    eqn1_file = fopen(str1, 'wt');
    fprintf(eqn1_file, '%s\n', char(eqns_(1)));
    fclose(eqn1_file);
    
    str2 = sprintf('eqn2_singularität');
    eqn2_file = fopen(str2, 'wt');
    fprintf(eqn2_file, '%s\n', char(eqns_(2)));
    fclose(eqn2_file);
    
    str3 = sprintf('eqn3_singularität');
    eqn3_file = fopen(str3, 'wt');
    fprintf(eqn3_file, '%s\n', char(eqns_(3)));
    fclose(eqn3_file);
