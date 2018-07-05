function [bending_stiffness_derivative] = GetBendingStiffnessDeriv()
    syms C q_0 A k D_damp q_1_dot q_2_dot q_3_dot v_1_dot v_2_dot v_3_dot phi theta kappa r_b q_1 q_2 q_3 p_1 p_2 p_3 k_spline;
    
    d_theta_dq_q1 = -(q_2 - 2*q_1 + q_3)/(3*r_b*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));

    d_theta_dq_q2 = -(q_1 - 2*q_2 + q_3)/(3*r_b*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));
    d_theta_dq_q3 = -(q_1 + q_2 - 2*q_3)/(3*r_b*(q_1^2 - q_1*q_2 - q_1*q_3 + q_2^2 - q_2*q_3 + q_3^2)^(1/2));

    d_M_1_bend_dq1 = C*d_theta_dq_q1;
    d_M_1_bend_dq2 = C*d_theta_dq_q2;
    d_M_1_bend_dq3 = C*d_theta_dq_q3;

    bending_stiffness_derivative(1,1) = d_M_1_bend_dq1;
    bending_stiffness_derivative(2,1) = d_M_1_bend_dq2;
    bending_stiffness_derivative(3,1) = d_M_1_bend_dq3;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
end

