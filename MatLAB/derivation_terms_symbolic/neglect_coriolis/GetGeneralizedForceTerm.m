
function [generalized_force_term] = GetGeneralizedForceTerm(w_H_h, dw_H_h__dq1, dw_H_h__dq2, dw_H_h__dq3)
    
    syms m q_0 A k D_damp_spline D_damp q_1_dot q_2_dot q_3_dot v_1_dot v_2_dot v_3_dot phi theta kappa r_b q_1 q_2 q_3 p_1 p_2 p_3 k_spline;
    
    l_bar = theta/kappa;
    l_bar_dot = (q_1_dot + q_2_dot + q_3_dot)/3;
    
    % Define forces on muscles in head frame
    h_F_active_muscle_1 = p_1*A*[0; 0; 1; 0];
    h_F_active_muscle_2 = p_2*A*[0; 0; 1; 0];
    h_F_active_muscle_3 = p_3*A*[0; 0; 1; 0];

    h_F_passive_muscle_1 = - (k*(q_1 - q_0) + D_damp*q_1_dot)*[0; 0; 1; 0];  
    h_F_passive_muscle_2 = - (k*(q_2 - q_0) + D_damp*q_2_dot)*[0; 0; 1; 0];
    h_F_passive_muscle_3 = - (k*(q_3 - q_0) + D_damp*q_3_dot)*[0; 0; 1; 0];

    h_F_muscle_1 = h_F_active_muscle_1 + h_F_passive_muscle_1; 
    h_F_muscle_2 = h_F_active_muscle_2 + h_F_passive_muscle_2;
    h_F_muscle_3 = h_F_active_muscle_3 + h_F_passive_muscle_3;
    
    % Define spline spring force.
    h_F_spline = - (k_spline*(l_bar - q_0) + D_damp_spline*l_bar_dot)*[0; 0; 1; 0];
    
    % Transform all forces w.r.t to world frame.
    w_F_muscle_1 = w_H_h*h_F_muscle_1;
    w_F_muscle_2 = w_H_h*h_F_muscle_2;
    w_F_muscle_3 = w_H_h*h_F_muscle_3;
    w_F_spline = w_H_h*h_F_spline;
    
    % Define position of muscle in head frame
    h_r_muscle_1 = [r_b; 0; 0; 1];
    h_r_muscle_2 = [- r_b/2; (sqrt(3)/2)*r_b; 0; 1]; 
    h_r_muscle_3 = [- r_b/2; - (sqrt(3)/2)*r_b; 0; 1]; 
    h_r_m = [0; 0; 0; 1];
    
    term_q1_chi_1 = (w_F_muscle_1.')*dw_H_h__dq1*h_r_muscle_1; 
    term_q1_chi_2 = (w_F_muscle_2.')*dw_H_h__dq1*h_r_muscle_2;
    term_q1_chi_3 = (w_F_muscle_3.')*dw_H_h__dq1*h_r_muscle_3;
    term_q1_spline = (w_F_spline.')*dw_H_h__dq1*h_r_m;
    term_q1 = term_q1_spline + term_q1_chi_1 + term_q1_chi_2 + term_q1_chi_3; 

    term_q2_chi_1 = w_F_muscle_1.'*dw_H_h__dq2*h_r_muscle_1;
    term_q2_chi_2 = w_F_muscle_2.'*dw_H_h__dq2*h_r_muscle_2;
    term_q2_chi_3 = w_F_muscle_3.'*dw_H_h__dq2*h_r_muscle_3;
    term_q2_spline = (w_F_spline.')*dw_H_h__dq2*h_r_m;
    term_q2 = term_q2_spline + term_q2_chi_1 + term_q2_chi_2 + term_q2_chi_3;

    term_q3_chi_1 = w_F_muscle_1.'*dw_H_h__dq3*h_r_muscle_1;
    term_q3_chi_2 = w_F_muscle_2.'*dw_H_h__dq3*h_r_muscle_2;
    term_q3_chi_3 = w_F_muscle_3.'*dw_H_h__dq3*h_r_muscle_3;
    term_q3_spline = (w_F_spline.')*dw_H_h__dq3*h_r_m;
    term_q3 = term_q3_spline + term_q3_chi_1 + term_q3_chi_2 + term_q3_chi_3;

    generalized_force_term = [term_q1; term_q2; term_q3];
    
end 