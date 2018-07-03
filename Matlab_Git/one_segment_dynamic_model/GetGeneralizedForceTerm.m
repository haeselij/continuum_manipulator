
function [generalized_force_term] = GetGeneralizedForceTerm(w_H_1h, dw_H_1h__dq1,dw_H_1h__dq2, dw_H_1h__dq3, q_1, q_2, q_3, p_1, p_2, p_3, q_dot)

    r_b = 0.026;
    A = 0.0183*0.0183*pi;  
    k = 20000;% "real value" -> 1108.63; %[N/m]
    q_0 = 0.3014;
    D_damp = 100;
    roundn = @(x,n) round(x*10^n)./10^n;

    % Define position of muscle in head frame
    h1_r_B1 = [r_b; 0; 0; 1];
    h1_r_B2 = [-(r_b)/2; (sqrt(3)/2)*(r_b); 0; 1]; 
    h1_r_B3 =[-(r_b)/2; -(sqrt(3)/2)*(r_b); 0; 1]; 

    % Define forces on muscles in head frame
    h1_F_active_B1 = p_1*A*[0; 0; 1; 1];
    h1_F_active_B2 = p_2*A*[0; 0; 1; 1];
    h1_F_active_B3 = p_3*A*[0; 0; 1; 1];

    h1_F_passive_B1 = -(k*(q_1 - q_0) + D_damp*q_dot(1))*[0; 0; 1; 1];  
    h1_F_passive_B2 = -(k*(q_2 - q_0) + D_damp*q_dot(2))*[0; 0; 1; 1];
    h1_F_passive_B3 = -(k*(q_3 - q_0) + D_damp*q_dot(3))*[0; 0; 1; 1];

    h1_F_B1 = h1_F_active_B1 + h1_F_passive_B1; 
    h1_F_B2 = h1_F_active_B2 + h1_F_passive_B2;
    h1_F_B3 = h1_F_active_B3 + h1_F_passive_B3;

    % Transform forces w.r.t world frame
    w_F_B1 = w_H_1h*h1_F_B1;
    w_F_B2 = w_H_1h*h1_F_B2;
    w_F_B3 = w_H_1h*h1_F_B3;

    term_q1_chi_1 = (w_F_B1.')*dw_H_1h__dq1*h1_r_B1; 
    term_q1_chi_2 = (w_F_B2.')*dw_H_1h__dq1*h1_r_B2;
    term_q1_chi_3 = (w_F_B3.')*dw_H_1h__dq1*h1_r_B3;
    term_q1 = term_q1_chi_1 + term_q1_chi_2 + term_q1_chi_3; 
    term_q1 = roundn(term_q1,10);


    term_q2_chi_1 = w_F_B1.'*dw_H_1h__dq2*h1_r_B1;
    term_q2_chi_2 = w_F_B2.'*dw_H_1h__dq2*h1_r_B2;
    term_q2_chi_3 = w_F_B3.'*dw_H_1h__dq2*h1_r_B3;
    term_q2 = term_q2_chi_1 + term_q2_chi_2 + term_q2_chi_3;
    term_q2 = roundn(term_q2,10);

    term_q3_chi_1 = w_F_B1.'*dw_H_1h__dq3*h1_r_B1;
    term_q3_chi_2 = w_F_B2.'*dw_H_1h__dq3*h1_r_B2;
    term_q3_chi_3 = w_F_B3.'*dw_H_1h__dq3*h1_r_B3;
    term_q3 = term_q3_chi_1 + term_q3_chi_2 + term_q3_chi_3;
    term_q3 = roundn(term_q3,10);

    generalized_force_term = [term_q1; term_q2; term_q3];
end 
