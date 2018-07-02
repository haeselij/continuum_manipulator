syms q_1 q_2 q_3 r_b
assume(q_1,'real');
assume(q_2,'real');
assume(q_3,'real');
assume(v_1,'real');
assume(v_2,'real');
assume(v_3,'real');

theta = 2/3*(sqrt(q_1^2 + q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3 - q_2*q_3))/r_b;
l_bar = (q_1 + q_2 + q_3)/3;
phi = atan2((sqrt(3)*(q_3 - q_2)),(q_2 + q_3 - 2*q_1));
kappa = theta/l_bar;
w_H_1h_1 = [cos(phi)^2*(cos(theta) - 1) + 1, sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)*sin(theta), (cos(phi)*(1 - cos(theta)))/kappa];
w_H_1h_2 = [sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)^2*(1 - cos(theta)) + cos(theta), sin(phi)*sin(theta), (sin(phi)*(1 - cos(theta)))/kappa];
w_H_1h_3 = [-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta), sin(theta)/kappa];
w_H_1h_4 = [0, 0, 0, 1];
w_H_1h = [w_H_1h_1; w_H_1h_2; w_H_1h_3; w_H_1h_4];

d_w_H_1h__dq1 =  diff(w_H_1h,q_1);
d_w_H_1h__dq2 =  diff(w_H_1h,q_2);
d_w_H_1h__dq3 =  diff(w_H_1h,q_3);

k_sys = [phi; kappa; l_bar];
q_ =[q_1 q_2 q_3];
for i=1:3
d_phi_dq(i) = diff(phi,q_(i));
d_kappa_dq(i) = diff(kappa,q_(i));
d_l_bar_dq(i) = diff(l_bar,q_(i));
d2_phi_dq1q(i) = diff(phi,q_(i));
for k=1:3
d2_phi_ddq(i,k) = diff(d_phi_dq(i),q_(k));
d2_kappa_ddq(i,k) = diff(d_kappa_dq(i),q_(k));
d2_l_bar_ddq(i,k) = diff(d_l_bar_dq(i),q_(k));
end
end

d_k_dq1 = diff(k_sys,q_1);
d_k_dq2 = diff(k_sys,q_2);
d_k_dq3 = diff(k_sys,q_3);

d2_k_dq1dq1 =  diff(d_k_dq1,q_1);
d2_k_dq1dq2 = diff(d_k_dq1,q_2);
d2_k_dq1dq3 = diff(d_k_dq1,q_3);
d2_k_dq2dq2 = diff(d_k_dq2,q_2);
d2_k_dq2dq3 = diff(d_k_dq2,q_3);
d2_k_dq3dq3 = diff(d_k_dq3,q_3);


