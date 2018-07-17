% Main script in order to derive Lagrange equations symbolically. Corresponding function
% files are SetMatrices, GetBendingStiffnessDeriv, GetGeneralizedForceTerm,
% GetKineticEnergyTerm and LagrangeEquation.

%%
% Main section, where all terms of Lagrange equarions are collected and
% saved. The results are the symbolic equations.

    % !!!!!! Equations are derived for two cases, one is the general case and the other
    % one is the case where singularity occurs. Depednig on for which case
    % the equations have to be derived 'calculate_singularity_case' is set
    % HERE to true for singularity case and false for general one. !!!!!!
    calculate_singularity_case = true;

    [w_H_h_, w_H_2dot_, dw_H_h__dq1_, dw_H_h__dq2_, dw_H_h__dq3_]= SetMatrices(calculate_singularity_case);
    [bending_stiffness_derivative_] = GetBendingStiffnessDeriv();

    [generalized_force_term_] = GetGeneralizedForceTerm(w_H_h_, dw_H_h__dq1_, dw_H_h__dq2_, dw_H_h__dq3_);
    [kinetic_energy_term_] = GetKineticEnergyTerm(dw_H_h__dq1_, dw_H_h__dq2_, dw_H_h__dq3_, w_H_2dot_);
    [potential_energy_term_] = GetPotentialEnergyTerm(dw_H_h__dq1_, dw_H_h__dq2_, dw_H_h__dq3_, bending_stiffness_derivative_);

    [eqns_] = LagrangeEquation(generalized_force_term_, kinetic_energy_term_, potential_energy_term_);

%%
% Print section, where equations are printed in text files. Depending one which one we
% are calculating, the filename (defined under variable 'strX') has to
% be adapted.

    if(calculate_singularity_case)

        str1 = sprintf('eqn1_singularity');
        eqn1_file = fopen(str1, 'wt');
        fprintf(eqn1_file, '%s\n', char(eqns_(1)));
        fclose(eqn1_file);

        str2 = sprintf('eqn2_singularity');
        eqn2_file = fopen(str2, 'wt');
        fprintf(eqn2_file, '%s\n', char(eqns_(2)));
        fclose(eqn2_file);

        str3 = sprintf('eqn3_singularity');
        eqn3_file = fopen(str3, 'wt');
        fprintf(eqn3_file, '%s\n', char(eqns_(3)));
        fclose(eqn3_file);

    else

        str1 = sprintf('eqn1');
        eqn1_file = fopen(str1, 'wt');
        fprintf(eqn1_file, '%s\n', char(eqns_(1)));
        fclose(eqn1_file);

        str2 = sprintf('eqn2');
        eqn2_file = fopen(str2, 'wt');
        fprintf(eqn2_file, '%s\n', char(eqns_(2)));
        fclose(eqn2_file);

        str3 = sprintf('eqn3');
        eqn3_file = fopen(str3, 'wt');
        fprintf(eqn3_file, '%s\n', char(eqns_(3)));
        fclose(eqn3_file);

    end
