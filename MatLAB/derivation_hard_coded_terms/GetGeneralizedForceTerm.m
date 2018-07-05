
function [generalized_force_term] = GetGeneralizedForceTerm(w_H_1h, dw_H_1h__dq1, dw_H_1h__dq2, dw_H_1h__dq3)
    syms m q_0 A k D_damp_spline D_damp q_1_dot q_2_dot q_3_dot v_1_dot v_2_dot v_3_dot phi theta kappa r_b q_1 q_2 q_3 p_1 p_2 p_3 k_spline;
    
    l_bar = theta/kappa;
    l_bar_dot = (q_1_dot + q_2_dot + q_3_dot)/3;
    
    % Define position of muscle in head frame
    h1_r_B1 = [r_b; 0; 0; 1];
    h1_r_B2 = [-(r_b)/2; (sqrt(3)/2)*(r_b); 0; 1]; 
    h1_r_B3 = [-(r_b)/2; -(sqrt(3)/2)*(r_b); 0; 1]; 
    h1_r_m = [0; 0; 0; 1];

    % Define forces on muscles in head frame
    h1_F_active_B1 = p_1*A*[0; 0; 1; 0];
    h1_F_active_B2 = p_2*A*[0; 0; 1; 0];
    h1_F_active_B3 = p_3*A*[0; 0; 1; 0];

    h1_F_passive_B1 = -(k*(q_1 - q_0) + D_damp*q_1_dot)*[0; 0; 1; 0];  
    h1_F_passive_B2 = -(k*(q_2 - q_0) + D_damp*q_2_dot)*[0; 0; 1; 0];
    h1_F_passive_B3 = -(k*(q_3 - q_0) + D_damp*q_3_dot)*[0; 0; 1; 0];

    h1_F_B1 = h1_F_active_B1 + h1_F_passive_B1; 
    h1_F_B2 = h1_F_active_B2 + h1_F_passive_B2;
    h1_F_B3 = h1_F_active_B3 + h1_F_passive_B3;
    
    h1_F_spline = -(k_spline*(l_bar - q_0) + D_damp_spline*l_bar_dot)*[0; 0; 1; 0];
    
    % Transform forces w.r.t world frame
    w_F_B1 = w_H_1h*h1_F_B1;
    w_F_B2 = w_H_1h*h1_F_B2;
    w_F_B3 = w_H_1h*h1_F_B3;
    w_F_spline = w_H_1h*h1_F_spline;
    
    term_q1_chi_1 = (w_F_B1.')*dw_H_1h__dq1*h1_r_B1; 
    term_q1_chi_2 = (w_F_B2.')*dw_H_1h__dq1*h1_r_B2;
    term_q1_chi_3 = (w_F_B3.')*dw_H_1h__dq1*h1_r_B3;
    term_q1_spline = (w_F_spline.')*dw_H_1h__dq1*h1_r_m;
    term_q1 = term_q1_spline + term_q1_chi_1 + term_q1_chi_2 + term_q1_chi_3; 

    term_q2_chi_1 = w_F_B1.'*dw_H_1h__dq2*h1_r_B1;
    term_q2_chi_2 = w_F_B2.'*dw_H_1h__dq2*h1_r_B2;
    term_q2_chi_3 = w_F_B3.'*dw_H_1h__dq2*h1_r_B3;
    term_q2_spline = (w_F_spline.')*dw_H_1h__dq2*h1_r_m;
    term_q2 = term_q2_spline + term_q2_chi_1 + term_q2_chi_2 + term_q2_chi_3;

    term_q3_chi_1 = w_F_B1.'*dw_H_1h__dq3*h1_r_B1;
    term_q3_chi_2 = w_F_B2.'*dw_H_1h__dq3*h1_r_B2;
    term_q3_chi_3 = w_F_B3.'*dw_H_1h__dq3*h1_r_B3;
    term_q3_spline = (w_F_spline.')*dw_H_1h__dq3*h1_r_m;
    term_q3 = term_q3_spline + term_q3_chi_1 + term_q3_chi_2 + term_q3_chi_3;

    generalized_force_term = [term_q1; term_q2; term_q3];
end 