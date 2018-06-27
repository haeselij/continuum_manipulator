dt = 0.01;
q_1 =0.31;
q_2 =0.3;
q_3 = 0.32;
v_1 = 0.2;
v_2 = 0;
v_3 = -0.01;
r_s = 0.01;
p_1 =  0.1;
p_2 = p_1;
p_3 = p_2;
r_b = 0.035;
A = 0.0183*0.0183*pi;  
k=  0.0111;
q_0 = 0.3014;
% Define position of muscle in head frame
h1_r_B1 = [0; r_b/2; 0; 0];
h1_r_B2 = [(sqrt(3)/2)*(r_b/2); -(r_b/2)/2; 0; 0];
h1_r_B3 = [-(sqrt(3)/2)*(r_b/2); -(r_b/2)/2; 0; 0];

% Define forces on muscles in head frame
h1_F_active_B1 = p_1*A*[0; 0; 1; 0];
h1_F_active_B2 = p_2*A*[0; 0; 1; 0];
h1_F_active_B3 = p_3*A*[0; 0; 1; 0];

h1_F_passive_B1 = -k*(q_1 - q_0)*[0; 0; 1; 0];  
h1_F_passive_B2 = -k*(q_2 - q_0)*[0; 0; 1; 0];
h1_F_passive_B3 = -k*(q_3 - q_0)*[0; 0; 1; 0];

h1_F_B1 = h1_F_active_B1 + h1_F_passive_B1; 
h1_F_B2 = h1_F_active_B2 + h1_F_passive_B2;
h1_F_B3 = h1_F_active_B3 + h1_F_passive_B3;

% Transform forces w.r.t world frame
w_F_B1 = w_H_1h*h1_F_B1;
w_F_B2 = w_H_1h*h1_F_B2;
w_F_B3 = w_H_1h*h1_F_B3;

term_q1_chi_1 = w_F_B1*dw_H_1h__dq1*h1_r_B1;
term_q1_chi_2 = w_F_B2*dw_H_1h__dq1*h1_r_B2;
term_q1_chi_3 = w_F_B3*dw_H_1h__dq1*h1_r_B3;
term_q1 = term_q1_chi_1 + term_q1_chi_2 + term_q1_chi_3;

term_q2_chi_1 = w_F_B1*dw_H_1h__dq2*h1_r_B1;
term_q2_chi_2 = w_F_B2*dw_H_1h__dq2*h1_r_B2;
term_q2_chi_3 = w_F_B3*dw_H_1h__dq2*h1_r_B3;
term_q2 = term_q2_chi_1 + term_q2_chi_2 + term_q2_chi_3;


term_q3_chi_1 = w_F_B1*dw_H_1h__dq3*h1_r_B1;
term_q3_chi_2 = w_F_B2*dw_H_1h__dq3*h1_r_B2;
term_q3_chi_3 = w_F_B3*dw_H_1h__dq3*h1_r_B3;
term_q3 = term_q2_chi_1 + term_q3_chi_2 + term_q3_chi_3;

generalized_force_term = [term_q1; term_q2; term_q3];